import rebound as rb
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
import pymcel as pym
import plotly.graph_objects as go
import plotly.express as pxi
from matplotlib.animation import FuncAnimation

from astropy.time import Time
import astropy.units as u 
from datetime import datetime, timedelta


# Funcion para integrar y dejar los vectores conrespecto a la posicion de la tierra

def integracion_Apo_tierr(sim, ts):
  sim_int = sim.copy()  # Creamos una copia de la simulación original para no modificarla durante el proceso
  
  n_cuerpos = len(sim.particles)
  rs = np.zeros((n_cuerpos, len(ts), 3))
  vs = np.zeros((n_cuerpos, len(ts), 3))
  rs_Tier = np.zeros((n_cuerpos, len(ts), 3))
  vs_Tier = np.zeros((n_cuerpos, len(ts), 3))
  dE = np.zeros(len(ts))

  E0 = sim_int.energy()
  #Integración
  for i, t in enumerate(ts):
      sim_int.integrate(t)

      # Guardar energía
      dE[i] = abs(sim_int.energy() - E0)/E0

      pos_jupiter = sim_int.particles[0].xyz
      vel_jupiter = sim_int.particles[0].vxyz

      for j in range(n_cuerpos):
        # Acceso directo a los datos de la partícula
        p = sim_int.particles[j]
        rs[j, i] = p.xyz
        vs[j, i] = p.vxyz
        rs_Tier[j, i] = rs[j, i] - rs[0, i]
        vs_Tier[j, i] = vs[j, i] - vs[0, i]

  return rs, vs, rs_Tier, vs_Tier, dE



# rutina para calcular la distancia minima de Apophis a la tierra y el tiempo

def distancia_minima(sim, ts, fecha_inicial, names):
  
  sim0 = sim.copy()  # Creamos una copia de la simulación original para no modificarla durante el proceso

  dis = []
  ts_min = []
  ts_d = []

  ts = np.asarray(ts, dtype=float)
  if ts.ndim != 1 or len(ts) < 2:
    raise ValueError("ts debe ser un arreglo 1D con al menos dos instantes de tiempo.")

  # Refinamos el minimo en una ventana local para reducir errores por muestreo grueso.
  ventana_refinamiento = max((ts[-1] - ts[0]) / 200.0, 0.01)
  n_refinado = max(5000, len(ts) * 5)

  for i in range(len(names)):
    sim1 = sim0.copy()

    j = 0


    while j <= i:
      sim1.add(names[j], date=fecha_inicial)
      j += 1
    
    n_bodies = len(sim1.particles)

    print(i, names[i], n_bodies)

    rs, vs, rs_Tier, vs_Tier, dE = integracion_Apo_tierr(sim1, ts)

    # Calculamos la distancia mínima entre Apophis y la Tierra
    distancias = np.linalg.norm(rs_Tier[1] - rs_Tier[0], axis=1)
    i_min = np.argmin(distancias)
    t_cercano = ts[i_min]

    t_ini = max(ts[0], t_cercano - ventana_refinamiento)
    t_fin = min(ts[-1], t_cercano + ventana_refinamiento)
    ts_fino = np.linspace(t_ini, t_fin, n_refinado)

    if len(ts_fino) >= 2 and t_fin > t_ini:
      rs_fino, vs_fino, rs_Tier_fino, vs_Tier_fino, dE_fino = integracion_Apo_tierr(sim1, ts_fino)
      distancias_finas = np.linalg.norm(rs_Tier_fino[1] - rs_Tier_fino[0], axis=1)
      i_min_fino = np.argmin(distancias_finas)
      S = distancias_finas[i_min_fino]
      tiempo_minimo = ts_fino[i_min_fino]
    else:
      S = distancias[i_min]
      tiempo_minimo = t_cercano
    



    distancia_minimakm = S * 149597870.7
    tdias = tiempo_minimo * 365.25
    print(f"Distancia mínima entre Apophis y la Tierra: {distancia_minimakm:.2f} km")
    print(tdias)
    t0 = Time(fecha_inicial, scale="utc")
    t1 = t0 + tdias * u.day

    dis.append(distancia_minimakm)
    ts_min.append(t1)
    ts_d.append(tdias)

  return rs, vs, dis, ts_min, ts_d


def plot_trayectoria_apophis(rs):

  color = ['red', 'yellow']

  plt.figure(figsize=(10, 8))
  plt.plot(rs[:, 0], rs[:, 1], label='Apophis', color=f'{color[0]}')
  plt.scatter(rs[0, 0], rs[0, 1], color=f'{color[0]}', marker='o', s=100, label='Apophis (inicio)')
  plt.scatter(rs[-1, 0], rs[-1, 1], color=f'{color[0]}', marker='x', s=100, label='Apophis (final)')
  plt.scatter(0, 0, color='blue', marker='o', s=100, label='Tierra (inicio)')
  plt.xlabel('X (AU)')
  plt.ylabel('Y (AU)')
  plt.title('Trayectorias de los cuerpos en el sistema solar respecto a la Tierra')
  plt.legend()
  plt.axis('equal')
  plt.grid()
  plt.show()


def construir_etiquetas_tiempo(fecha_inicial, ts, utc_offset_horas=-5):
  """Construye etiquetas de tiempo UTC y UTC offset para cada frame.

  Parametros
  ----------
  fecha_inicial : str o astropy.time.Time
      Fecha inicial de la simulacion.
  ts : array-like
      Arreglo de tiempos en anos.
  utc_offset_horas : float, opcional
      Offset horario respecto a UTC para la segunda etiqueta.

  Retorna
  -------
  list[tuple[str, str]]
      Lista con pares (texto_utc, texto_utc_offset) por frame.
  """
  ts = np.asarray(ts, dtype=float)
  t0 = fecha_inicial if isinstance(fecha_inicial, Time) else Time(fecha_inicial, scale="utc")
  dias = ts * 365.25 * u.day
  tiempos_utc = t0 + dias
  tiempos_local = tiempos_utc + (utc_offset_horas * u.hour)

  etiquetas = []
  for t_utc, t_loc in zip(tiempos_utc, tiempos_local):
    etiquetas.append(
      (
        f"UTC: {t_utc.to_datetime().strftime('%Y-%m-%d %H:%M:%S')}",
        f"UTC{utc_offset_horas:+g}: {t_loc.to_datetime().strftime('%Y-%m-%d %H:%M:%S')}",
      )
    )
  return etiquetas


def animacion_2d_tierra_apophis(
  rs,
  intervalo_ms=40,
  trail=250,
  show=True,
  fecha_inicial=None,
  ts=None,
  utc_offset_horas=-5,
  max_frames=900,
):
  """Anima trayectorias 2D baricentricas y relativas Tierra-Apophis.

  Parametros
  ----------
  rs : np.ndarray
      Arreglo de posiciones con forma (2, ts, 3):
      rs[0] = Tierra respecto al baricentro,
      rs[1] = Apophis respecto al baricentro.
  intervalo_ms : int, opcional
      Intervalo entre frames en milisegundos.
  trail : int o None, opcional
      Longitud del rastro mostrado. Si es None, dibuja toda la trayectoria
      acumulada hasta el frame actual.
  show : bool, opcional
      Si es True, muestra la figura al final.
  fecha_inicial : str o astropy.time.Time, opcional
      Fecha inicial para etiquetar el tiempo en la animacion.
  ts : np.ndarray, opcional
      Arreglo de tiempos en anos asociado a cada frame. Debe tener largo n_t.
  utc_offset_horas : float, opcional
      Offset horario adicional a mostrar (por defecto UTC-5).
    max_frames : int, opcional
      Maximo de frames que se usaran en la animacion. Si la serie temporal
      es mas larga, se muestrea de forma uniforme para acelerar la generacion.

  Retorna
  -------
  fig : matplotlib.figure.Figure
  ani : matplotlib.animation.FuncAnimation
  """
  rs = np.asarray(rs)

  if rs.ndim != 3 or rs.shape[0] != 2 or rs.shape[2] != 3:
    raise ValueError(
      "rs debe tener forma (2, ts, 3): [Tierra, Apophis] x tiempo x (x,y,z)."
    )

  n_t = rs.shape[1]
  if n_t < 2:
    raise ValueError("rs debe contener al menos 2 instantes de tiempo.")

  etiquetas_tiempo = None
  if fecha_inicial is not None or ts is not None:
    if fecha_inicial is None or ts is None:
      raise ValueError("Debes proporcionar ambos parametros: fecha_inicial y ts.")
    ts = np.asarray(ts)
    if ts.ndim != 1:
      raise ValueError("ts debe ser un arreglo 1D de tiempos en anos.")
    if len(ts) != n_t:
      raise ValueError("ts debe tener el mismo numero de frames que rs.shape[1].")
    etiquetas_tiempo = construir_etiquetas_tiempo(fecha_inicial, ts, utc_offset_horas)

  if max_frames is None or max_frames <= 0:
    frame_indices = np.arange(n_t)
  elif n_t <= max_frames:
    frame_indices = np.arange(n_t)
  else:
    frame_indices = np.unique(np.linspace(0, n_t - 1, int(max_frames), dtype=int))

  tierra_xy = rs[0, :, :2]
  apophis_xy = rs[1, :, :2]
  apo_rel_xy = (rs[1] - rs[0])[:, :2]

  fig, axs = plt.subplots(1, 2, figsize=(14, 6))

  # Panel izquierdo: marco baricentrico
  ax0 = axs[0]
  ax0.set_title("Baricentro del Sistema Solar")
  ax0.set_xlabel("X (AU)")
  ax0.set_ylabel("Y (AU)")
  ax0.grid(alpha=0.3)
  ax0.scatter(0, 0, c="gold", s=120, marker="*", label="Sol (0,0)")

  line_tierra, = ax0.plot([], [], lw=1.8, c="tab:blue", label="Tierra")
  line_apophis, = ax0.plot([], [], lw=1.8, c="tab:red", label="Apophis")
  p_tierra, = ax0.plot([], [], marker="o", c="tab:blue", ms=6)
  p_apophis, = ax0.plot([], [], marker="o", c="tab:red", ms=6)
  ax0.legend(loc="best")

  # Panel derecho: Apophis respecto a Tierra
  ax1 = axs[1]
  ax1.set_title("Apophis respecto a la Tierra")
  ax1.set_xlabel("X (AU)")
  ax1.set_ylabel("Y (AU)")
  ax1.grid(alpha=0.3)
  ax1.scatter(0, 0, c="tab:blue", s=40, marker="o", label="Tierra (origen)")

  line_rel, = ax1.plot([], [], lw=1.8, c="tab:orange", label="Apophis - Tierra")
  p_rel, = ax1.plot([], [], marker="o", c="tab:orange", ms=6)
  ax1.legend(loc="best")

  time_text = None
  if etiquetas_tiempo is not None:
    time_text = ax1.text(
      0.02,
      0.02,
      "",
      transform=ax1.transAxes,
      ha="left",
      va="bottom",
      fontsize=10,
      bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none"},
    )

  # Limites con margen para que ambos paneles queden legibles.
  xy_all_0 = np.vstack([tierra_xy, apophis_xy, np.array([[0.0, 0.0]])])
  x0_min, y0_min = xy_all_0.min(axis=0)
  x0_max, y0_max = xy_all_0.max(axis=0)
  dx0 = x0_max - x0_min
  dy0 = y0_max - y0_min
  pad0 = 0.08 * max(dx0, dy0, 1e-9)

  ax0.set_xlim(x0_min - pad0, x0_max + pad0)
  ax0.set_ylim(y0_min - pad0, y0_max + pad0)
  ax0.set_aspect("equal", adjustable="box")

  xy_all_1 = np.vstack([apo_rel_xy, np.array([[0.0, 0.0]])])
  x1_min, y1_min = xy_all_1.min(axis=0)
  x1_max, y1_max = xy_all_1.max(axis=0)
  dx1 = x1_max - x1_min
  dy1 = y1_max - y1_min
  pad1 = 0.08 * max(dx1, dy1, 1e-9)

  ax1.set_xlim(x1_min - pad1, x1_max + pad1)
  ax1.set_ylim(y1_min - pad1, y1_max + pad1)
  ax1.set_aspect("equal", adjustable="box")

  def _segmento(i):
    if trail is None:
      return 0, i + 1
    return max(0, i - int(trail) + 1), i + 1

  def init():
    line_tierra.set_data([], [])
    line_apophis.set_data([], [])
    p_tierra.set_data([], [])
    p_apophis.set_data([], [])
    line_rel.set_data([], [])
    p_rel.set_data([], [])
    artists = [line_tierra, line_apophis, p_tierra, p_apophis, line_rel, p_rel]
    if time_text is not None:
      txt_utc, txt_loc = etiquetas_tiempo[0]
      time_text.set_text(f"{txt_utc}    |    {txt_loc}")
      artists.append(time_text)
    return tuple(artists)

  def update(i):
    i0, i1 = _segmento(i)

    xt = tierra_xy[i0:i1, 0]
    yt = tierra_xy[i0:i1, 1]
    xa = apophis_xy[i0:i1, 0]
    ya = apophis_xy[i0:i1, 1]
    xr = apo_rel_xy[i0:i1, 0]
    yr = apo_rel_xy[i0:i1, 1]

    line_tierra.set_data(xt, yt)
    line_apophis.set_data(xa, ya)
    p_tierra.set_data([tierra_xy[i, 0]], [tierra_xy[i, 1]])
    p_apophis.set_data([apophis_xy[i, 0]], [apophis_xy[i, 1]])

    line_rel.set_data(xr, yr)
    p_rel.set_data([apo_rel_xy[i, 0]], [apo_rel_xy[i, 1]])

    artists = [line_tierra, line_apophis, p_tierra, p_apophis, line_rel, p_rel]
    if time_text is not None:
      txt_utc, txt_loc = etiquetas_tiempo[i]
      time_text.set_text(f"{txt_utc}    |    {txt_loc}")
      artists.append(time_text)

    return tuple(artists)

  ani = FuncAnimation(
    fig,
    update,
    frames=frame_indices,
    init_func=init,
    interval=intervalo_ms,
    blit=True,
    repeat=True,
  )

  plt.tight_layout()
  if show:
    plt.show()

  return fig, ani


def rotacion_peri_astronomico(r_rel, v_rel, Omega, I, omega):
  """Rota un vector de posición y velocidad desde el marco perifocal al inercial.

  Parametros
  ----------
  r_rel : array-like
      Vector de posición relativo (ts, 3).
  v_rel : array-like
      Vector de velocidad relativo (ts, 3).
  Omega : float
      Ascensión recta del nodo ascendente en radianes.
  I : float
      Inclinación orbital en radianes.
  omega : float
      Argumento del periastro en radianes.

  Retorna
  -------
  r_inercial : np.ndarray
      Vector de posición en el marco inercial (ts, 3).
  v_inercial : np.ndarray
      Vector de velocidad en el marco inercial (ts, 3).
  """
  # Convertir ángulos a radianes
  Omega_rad = Omega
  i_rad = I
  omega_rad = omega

  # Matriz de rotación perifocal a inercial
  R = np.array([
    [
      np.cos(Omega_rad) * np.cos(omega_rad) - np.sin(Omega_rad) * np.sin(omega_rad) * np.cos(i_rad),
      -np.cos(Omega_rad) * np.sin(omega_rad) - np.sin(Omega_rad) * np.cos(omega_rad) * np.cos(i_rad),
      np.sin(Omega_rad) * np.sin(i_rad)
    ],
    [
      np.sin(Omega_rad) * np.cos(omega_rad) + np.cos(Omega_rad) * np.sin(omega_rad) * np.cos(i_rad),
      -np.sin(Omega_rad) * np.sin(omega_rad) + np.cos(Omega_rad) * np.cos(omega_rad) * np.cos(i_rad),
      -np.cos(Omega_rad) * np.sin(i_rad)
    ],
    [
      np.sin(omega_rad) * np.sin(i_rad),
      np.cos(omega_rad) * np.sin(i_rad),
      np.cos(i_rad)
    ]
  ])

  r_inercial = R @ r_rel
  v_inercial = R @ v_rel

  return r_inercial, v_inercial
