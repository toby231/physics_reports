import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter
from mpl_toolkits.mplot3d import Axes3D
from model_functions import rk4

R = 47
C = 1e-6
scale_factor = C * R * 1e3
t0 = 0
tmax = 10 / scale_factor
h = 0.125e-3 / scale_factor  
N = int((tmax - t0) / h)
half = N // 2

# --- Run your model ---
Rv = 128.0
t, x, _, _, _ = rk4(Rv=Rv)
t = t[half:]
x = x[half:,:]

# --- Parameters ---
fps = 30                     # 🎥 30 FPS for smoother real-time feel
rotation_speed = 0.6         # degrees per frame (≈360° / 600 frames)
max_frames = 600             # total frames (~20 s)
step = max(1, len(t) // max_frames)
t = t[::step]
x = x[::step]

# --- Figure setup ---
fig = plt.figure(figsize=(7, 6))
ax = fig.add_subplot(111, projection='3d')

# Static axes (don’t warp)
ax.set_xlim(-1.2, 1.2)
ax.set_ylim(-1.2, 1.2)
ax.set_zlim(-0.6, 0.6)
ax.set_xlabel(r"$x$ (V)")
ax.set_ylabel(r"$V_1$ (V)")
ax.set_zlabel(r"$V_2$ (V)")
ax.set_title("3D Phase-Space Trajectory (30 FPS Loop)", fontsize=12)
ax.grid(True)

# Trajectory line + moving point
line, = ax.plot([], [], [], lw=0.5, color='royalblue')
point, = ax.plot([], [], [], 'o', color='orange', markersize=6)

# --- Init and update functions ---
def init():
    line.set_data([], [])
    line.set_3d_properties([])
    point.set_data([], [])
    point.set_3d_properties([])
    return line, point

def update(frame):
    line.set_data(x[:frame, 0], x[:frame, 1])
    line.set_3d_properties(x[:frame, 2])
    point.set_data(x[frame-1:frame, 0], x[frame-1:frame, 1])
    point.set_3d_properties(x[frame-1:frame, 2])
    # Rotate camera around z-axis
    ax.view_init(elev=25, azim=frame * rotation_speed)
    return line, point

# --- Build animation ---
ani = FuncAnimation(
    fig, update, frames=len(t),
    init_func=init, interval=1000/fps, blit=True
)

# --- Save MP4 ---
writer = FFMpegWriter(fps=fps, bitrate=2000)
ani.save(f"3D_phase_space_30fps_{Rv}.mp4", writer=writer, dpi=150)
print("✅ Saved 30 FPS MP4 → 3D_phase_space_30fps.mp4")

# --- Save looping GIF ---
ani.save(f"3D_phase_space_30fps_loop_{Rv}.gif", writer="pillow", fps=fps, dpi=120)
plt.close(fig)
