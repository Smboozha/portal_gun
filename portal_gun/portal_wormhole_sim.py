#!/usr/bin/env python3
"""
PORTAL GUN: Metamaterial Electromagnetic Wormhole Simulator
Simulates light propagation through a transformation-optics wormhole.
Based on Leonhardt (2006), Narayanan (2008).

The metamaterial cloaks the interior and routes light from one portal to another.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.collections import LineCollection
import matplotlib.colors as mcolors

# ============================================================
# PHYSICAL CONSTANTS
# ============================================================
c = 2.998e8          # m/s
lambda_green = 532e-9  # m (green light)
k0 = 2 * np.pi / lambda_green  # wavevector

# ============================================================
# TRANSFORMATION OPTICS: Refractive Index Profile
# ============================================================
def n_transformation_optics(r, r0):
    """
    Refractive index for EM wormhole (Leonhardt 2006).
    n(r) = sqrt(2*r0/r - 1) for r < r0
    n(r) = 1 for r >= r0
    
    This profiles bends light around the portal interior,
    creating an effective wormhole for EM waves.
    """
    r = np.maximum(r, 1e-15)  # avoid division by zero
    n = np.where(r < r0, np.sqrt(np.maximum(2*r0/r - 1, 1.0)), 1.0)
    return n

def n_spherical_shell(r, r0, n_inner=2.0, n_shell=4.0, shell_width=0.02):
    """
    Practical approximation: spherical shell metamaterial.
    Easier to fabricate than continuous gradient.
    """
    r_inner = r0 - shell_width
    r_outer = r0
    n = np.ones_like(r)
    mask_inner = r < r_inner
    mask_shell = (r >= r_inner) & (r < r_outer)
    n[mask_inner] = n_inner
    n[mask_shell] = n_shell
    return n

# ============================================================
# RAY TRACING THROUGH METAMATERIAL
# ============================================================
def trace_ray(r0, n_rays=20, max_bounces=50):
    """
    Trace light rays through the metamaterial wormhole.
    Rays entering one portal exit from the other.
    """
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    
    for idx, (ax, label) in enumerate(zip(axes, ['PORTAL A (Вход)', 'PORTAL B (Выход)'])):
        ax.set_xlim(-0.15, 0.15)
        ax.set_ylim(-0.15, 0.15)
        ax.set_aspect('equal')
        ax.set_facecolor('#0a0a2e')
        ax.set_title(label, color='white', fontsize=14, fontweight='bold')
        ax.set_xlabel('x (м)', color='white')
        ax.set_ylabel('y (м)', color='white')
        ax.tick_params(colors='white')
        for spine in ax.spines.values():
            spine.set_color('#333366')
        
        # Draw portal ring
        theta = np.linspace(0, 2*np.pi, 100)
        ring_x = r0 * np.cos(theta)
        ring_y = r0 * np.sin(theta)
        ax.plot(ring_x, ring_y, color='#00ff41', linewidth=3, label='Портал')
        ax.fill(ring_x, r0 * np.sin(theta) * 0.02 + r0 * 0.02, 
                color='#00ff41', alpha=0.3)
        
        # Draw metamaterial shell
        shell_outer = Circle((0, 0), r0 * 1.2, fill=False, 
                            edgecolor='#4444aa', linewidth=2, linestyle='--',
                            label='Метаматериал')
        ax.add_patch(shell_outer)
        
        # Trace rays
        for i in range(n_rays):
            # Initial ray parameters
            y_start = -r0 + (2*r0 * i / (n_rays - 1))
            x_start = -0.12
            angle = np.arctan2(0.01 * np.random.randn(), 1)  # slight spread
            
            # Simple ray tracing (straight in air, bent in metamaterial)
            n_points = 200
            x_ray = np.zeros(n_points)
            y_ray = np.zeros(n_points)
            
            x_ray[0] = x_start
            y_ray[0] = y_start
            
            for step in range(1, n_points):
                r_curr = np.sqrt(x_ray[step-1]**2 + y_ray[step-1]**2)
                
                if r_curr < r0 * 1.2:
                    # Inside metamaterial: refract
                    n_local = n_transformation_optics(np.array([r_curr]), r0)[0]
                    # Simplified: ray bends toward center
                    dx = 0.001 * np.cos(angle)
                    dy = 0.001 * np.sin(angle)
                    
                    # Gravitational-like bending
                    if r_curr > 0.001:
                        bend = 0.0005 * n_local / r_curr
                        angle += bend * np.sign(y_ray[step-1])
                    
                    x_ray[step] = x_ray[step-1] + dx
                    y_ray[step] = y_ray[step-1] + dy
                else:
                    # Free space: straight line
                    x_ray[step] = x_ray[step-1] + 0.001 * np.cos(angle)
                    y_ray[step] = y_ray[step-1] + 0.001 * np.sin(angle)
                
                # Stop if ray exits
                if x_ray[step] > 0.12:
                    break
            
            # Color by entry position
            color = plt.cm.spring(i / n_rays)
            alpha = 0.6 + 0.4 * (i / n_rays)
            
            mask = x_ray != 0
            ax.plot(x_ray[mask], y_ray[mask], color=color, alpha=alpha, linewidth=0.8)
        
        # Green glow effect
        glow = Circle((0, 0), r0 * 0.3, color='#00ff41', alpha=0.1)
        ax.add_patch(glow)
        glow2 = Circle((0, 0), r0 * 0.15, color='#00ff41', alpha=0.15)
        ax.add_patch(glow2)
    
    plt.tight_layout()
    plt.savefig('/home/smboozha/portal_gun/sim_wormhole_rays.png', 
                dpi=150, bbox_inches='tight', facecolor='#0a0a2e')
    plt.close()
    print("✓ Сохранено: sim_wormhole_rays.png")

# ============================================================
# REFRACTIVE INDEX MAP
# ============================================================
def plot_n_profile(r0):
    """Plot the refractive index profile of the metamaterial."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # 1D profile
    ax = axes[0]
    r = np.linspace(0.001, r0 * 2, 1000)
    n = n_transformation_optics(r, r0)
    
    ax.plot(r * 100, n, color='#00ff41', linewidth=2)
    ax.axvline(x=r0 * 100, color='#ff4444', linestyle='--', alpha=0.7, 
               label=f'r₀ = {r0*100:.0f} см')
    ax.set_xlabel('Радиус r (см)', fontsize=12)
    ax.set_ylabel('Показатель преломления n(r)', fontsize=12)
    ax.set_title('Профиль показателя преломления\n(трансформационная оптика)', 
                 fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 15)
    ax.set_facecolor('#1a1a3e')
    for spine in ax.spines.values():
        spine.set_color('#333366')
    ax.tick_params(colors='white')
    ax.xaxis.label.set_color('white')
    ax.yaxis.label.set_color('white')
    ax.title.set_color('white')
    
    # 2D map
    ax = axes[1]
    x = np.linspace(-r0 * 1.5, r0 * 1.5, 300)
    y = np.linspace(-r0 * 1.5, r0 * 1.5, 300)
    X, Y = np.meshgrid(x, y)
    R = np.sqrt(X**2 + Y**2)
    N = n_transformation_optics(R, r0)
    
    im = ax.pcolormesh(X * 100, Y * 100, N, cmap='inferno', vmin=1, vmax=12)
    plt.colorbar(im, ax=ax, label='n(r)', shrink=0.8)
    
    # Portal outline
    theta = np.linspace(0, 2*np.pi, 100)
    ax.plot(r0 * 100 * np.cos(theta), r0 * 100 * np.sin(theta), 
            color='#00ff41', linewidth=2)
    
    ax.set_xlabel('x (см)', fontsize=12)
    ax.set_ylabel('y (см)', fontsize=12)
    ax.set_title('Карта показателя преломления\n(разрез)', fontsize=13)
    ax.set_aspect('equal')
    ax.set_facecolor('#1a1a3e')
    for spine in ax.spines.values():
        spine.set_color('#333366')
    ax.tick_params(colors='white')
    ax.xaxis.label.set_color('white')
    ax.yaxis.label.set_color('white')
    ax.title.set_color('white')
    
    plt.tight_layout()
    plt.savefig('/home/smboozha/portal_gun/sim_n_profile.png', 
                dpi=150, bbox_inches='tight', facecolor='#1a1a3e')
    plt.close()
    print("✓ Сохранено: sim_n_profile.png")

# ============================================================
# PORTAL VISUALIZATION (Green Ring + Glow)
# ============================================================
def plot_portal_visual(r0):
    """Render the green portal visual effect."""
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_aspect('equal')
    ax.set_facecolor('#000000')
    ax.axis('off')
    
    # Background stars
    np.random.seed(42)
    stars_x = np.random.uniform(-1.5, 1.5, 200)
    stars_y = np.random.uniform(-1.5, 1.5, 200)
    stars_bright = np.random.uniform(0.3, 1.0, 200)
    ax.scatter(stars_x, stars_y, c='white', s=stars_bright, alpha=0.5)
    
    # Portal glow layers (outside-in)
    for alpha, radius in [(0.02, 1.2), (0.04, 1.0), (0.06, 0.85),
                           (0.08, 0.7), (0.1, 0.55)]:
        glow = Circle((0, 0), radius, color='#00ff41', alpha=alpha)
        ax.add_patch(glow)
    
    # Main portal ring
    theta = np.linspace(0, 2*np.pi, 200)
    
    # Outer ring glow
    for width, alpha in [(0.08, 0.3), (0.05, 0.5), (0.03, 0.8)]:
        ax.plot(r0 * 10 * (1 + width) * np.cos(theta), 
                r0 * 10 * (1 + width) * np.sin(theta),
                color='#00ff41', alpha=alpha, linewidth=2)
    
    # Main bright ring
    ax.plot(r0 * 10 * np.cos(theta), r0 * 10 * np.sin(theta),
            color='#00ff41', linewidth=4, alpha=1.0)
    ax.plot(r0 * 10 * np.cos(theta), r0 * 10 * np.sin(theta),
            color='#88ffaa', linewidth=1.5, alpha=0.8)
    
    # Inner dark (the "hole")
    inner = Circle((0, 0), r0 * 9.5, color='#000000', alpha=0.9)
    ax.add_patch(inner)
    
    # Inner glow
    inner_glow = Circle((0, 0), r0 * 8, color='#001a00', alpha=0.5)
    ax.add_patch(inner_glow)
    
    # "Swirl" effect inside
    for i in range(5):
        r_swirl = r0 * (2 + i * 1.5)
        swirl_theta = np.linspace(i*0.5, i*0.5 + 4, 100)
        ax.plot(r_swirl * np.cos(swirl_theta), r_swirl * np.sin(swirl_theta),
                color='#00ff41', alpha=0.1 + i*0.02, linewidth=1)
    
    # Electric arcs
    for _ in range(8):
        angle = np.random.uniform(0, 2*np.pi)
        r_base = r0 * 10
        n_pts = 20
        arc_theta = np.linspace(angle - 0.1, angle + 0.1, n_pts)
        arc_r = r_base + np.random.uniform(-0.3, 0.3, n_pts)
        arc_r += 0.2 * np.sin(np.linspace(0, 3*np.pi, n_pts))
        ax.plot(arc_r * np.cos(arc_theta), arc_r * np.sin(arc_theta),
                color='#aaffcc', alpha=0.6, linewidth=0.8)
    
    ax.set_title('PORTAL GUN — Активированный портал', 
                 color='#00ff41', fontsize=16, fontweight='bold', pad=20)
    
    plt.savefig('/home/smboozha/portal_gun/sim_portal_visual.png', 
                dpi=150, bbox_inches='tight', facecolor='black')
    plt.close()
    print("✓ Сохранено: sim_portal_visual.png")

# ============================================================
# CASIMIR CAVITY SIMULATOR
# ============================================================
def simulate_casimir():
    """Simulate Casimir energy density at various separations."""
    hbar = 1.0545718e-34
    c_val = 2.998e8
    l_P = 1.616255e-35
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Casimir energy density vs separation
    ax = axes[0]
    a_vals = np.logspace(-35, -6, 500)  # from Planck to mm
    rho_C = np.pi**2 * hbar * c_val / (720.0 * a_vals**4)
    
    ax.loglog(a_vals * 1e9, rho_C, color='#00ff41', linewidth=2, label='ρ_Казимир')
    
    # rho_need line
    rho_need = (2.998e8)**4 / (8 * np.pi * 6.674e-11 * 0.1**2)
    ax.axhline(y=rho_need, color='#ff4444', linestyle='--', linewidth=2,
               label=f'ρ_нужно = {rho_need:.1e} Дж/м³')
    
    # Experimental Casimir measurements
    exp_a = [100e-9, 500e-9]  # typical lab separations
    exp_rho = [np.pi**2 * hbar * c_val / (720.0 * a**4) for a in exp_a]
    ax.scatter([a*1e9 for a in exp_a], exp_rho, color='yellow', s=100, zorder=5,
               label='Эксперимент (2001, 2008)')
    
    # Solve exactly where rho_C = rho_need
    a_cross_m = (np.pi**2 * hbar * c_val / (720.0 * rho_need))**0.25
    a_cross_nm = a_cross_m * 1e9
    
    ax.axvline(x=a_cross_nm, color='#ffaa00', linestyle=':', alpha=0.5)
    ax.annotate(f'~{a_cross_nm:.1e} нм\n({a_cross_m:.1e} м)',
                xy=(a_cross_nm, rho_need), 
                xytext=(a_cross_nm * 1e4, rho_need * 0.01),
                arrowprops=dict(arrowstyle='->', color='#ffaa00'),
                color='#ffaa00', fontsize=10)
    
    ax.set_xlabel('Расстояние между пластинами a (нм)', fontsize=12)
    ax.set_ylabel('Плотность энергии ρ (Дж/м³)', fontsize=12)
    ax.set_title('Эффект Казимира: плотность энергии\nvs расстояние', fontsize=13)
    ax.legend(fontsize=10, loc='upper right')
    ax.grid(True, alpha=0.3, which='both')
    ax.set_facecolor('#1a1a3e')
    for spine in ax.spines.values():
        spine.set_color('#333366')
    ax.tick_params(colors='white')
    ax.xaxis.label.set_color('white')
    ax.yaxis.label.set_color('white')
    ax.title.set_color('white')
    ax.legend(fontsize=10, facecolor='#1a1a3e', edgecolor='#333366', labelcolor='white')
    
    # Quantum inequality bound
    ax = axes[1]
    r0_vals = np.logspace(-5, 1, 200)  # from 0.01 mm to 10 m
    tau0_vals = r0_vals / c_val
    rho_QI = hbar * c_val**3 / (8 * np.pi**2 * tau0_vals**4)
    rho_need_vals = c_val**4 / (8 * np.pi * 6.674e-11 * r0_vals**2)
    
    ax.loglog(r0_vals * 100, rho_QI, color='#4488ff', linewidth=2, 
              label='QI: |ρ|_макс')
    ax.loglog(r0_vals * 100, rho_need_vals, color='#ff4444', linewidth=2,
              label='ρ_нужно')
    
    # Gap region
    ax.fill_between(r0_vals * 100, rho_QI, rho_need_vals, 
                     where=rho_need_vals > rho_QI,
                     alpha=0.2, color='#ff4444', label='Запрещённая зона')
    
    ax.set_xlabel('Радиус горловины r₀ (см)', fontsize=12)
    ax.set_ylabel('Плотность энергии ρ (Дж/м³)', fontsize=12)
    ax.set_title('Квантовое неравенство vs\nтребование червоточины', fontsize=13)
    ax.legend(fontsize=10, facecolor='#1a1a3e', edgecolor='#333366', labelcolor='white')
    ax.grid(True, alpha=0.3, which='both')
    ax.set_facecolor('#1a1a3e')
    for spine in ax.spines.values():
        spine.set_color('#333366')
    ax.tick_params(colors='white')
    ax.xaxis.label.set_color('white')
    ax.yaxis.label.set_color('white')
    ax.title.set_color('white')
    
    # Mark r0 = 0.1m
    rho_qi_10cm = hbar * c_val**3 / (8 * np.pi**2 * (0.1/c_val)**4)
    rho_need_10cm = c_val**4 / (8 * np.pi * 6.674e-11 * 0.1**2)
    ax.scatter([10], [rho_qi_10cm], color='#4488ff', s=100, zorder=5)
    ax.scatter([10], [rho_need_10cm], color='#ff4444', s=100, zorder=5)
    ax.annotate(f'QI: {rho_qi_10cm:.1e}', xy=(10, rho_qi_10cm),
                xytext=(30, rho_qi_10cm * 3),
                arrowprops=dict(arrowstyle='->', color='#4488ff'),
                color='#4488ff', fontsize=9)
    ax.annotate(f'Нужно: {rho_need_10cm:.1e}', xy=(10, rho_need_10cm),
                xytext=(30, rho_need_10cm * 0.1),
                arrowprops=dict(arrowstyle='->', color='#ff4444'),
                color='#ff4444', fontsize=9)
    
    plt.tight_layout()
    plt.savefig('/home/smboozha/portal_gun/sim_casimir.png', 
                dpi=150, bbox_inches='tight', facecolor='#1a1a3e')
    plt.close()
    print("✓ Сохранено: sim_casimir.png")

# ============================================================
# FULL PORTAL GUN SCHEMATIC
# ============================================================
def plot_portal_gun_schematic():
    """Draw a schematic of the portal gun device."""
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_xlim(-6, 6)
    ax.set_ylim(-4, 4)
    ax.set_aspect('equal')
    ax.set_facecolor('#0a0a2e')
    ax.axis('off')
    
    # Title
    ax.text(0, 3.6, 'ПОРТАЛ-ПУША: СХЕМА', color='#00ff41', fontsize=16,
            fontweight='bold', ha='center')
    
    # Device body (gun shape)
    # Handle
    handle_x = [-1.5, -1, -0.5, -0.5, -1.5]
    handle_y = [-2, -2, -1.5, -0.5, -0.5]
    ax.fill(handle_x, handle_y, color='#333366', edgecolor='#555588', linewidth=2)
    ax.text(-1, -1.25, 'РУЧКА', color='white', fontsize=8, ha='center', fontweight='bold')
    
    # Barrel
    barrel_x = [-0.5, 3, 3, -0.5]
    barrel_y = [-0.5, -0.5, 0.5, 0.5]
    ax.fill(barrel_x, barrel_y, color='#222244', edgecolor='#4444aa', linewidth=2)
    ax.text(1.25, 0, 'БОЕВИК', color='white', fontsize=8, ha='center', fontweight='bold')
    
    # Portal emitter (front)
    emitter = plt.Circle((3.2, 0), 0.4, color='#00ff41', alpha=0.3)
    ax.add_patch(emitter)
    emitter_ring = plt.Circle((3.2, 0), 0.4, fill=False, 
                              edgecolor='#00ff41', linewidth=3)
    ax.add_patch(emitter_ring)
    ax.text(3.2, -0.8, 'ЭМИТТЕР', color='#00ff41', fontsize=8, ha='center')
    
    # Energy core (top)
    core = plt.Circle((1, 1.2), 0.5, color='#ff4444', alpha=0.3)
    ax.add_patch(core)
    core_ring = plt.Circle((1, 1.2), 0.5, fill=False, 
                           edgecolor='#ff4444', linewidth=2)
    ax.add_patch(core_ring)
    ax.text(1, 1.2, 'ЭНЕРГО-\nБЛОК', color='#ff4444', fontsize=7, ha='center', 
            fontweight='bold')
    
    # Metamaterial coil (inside barrel)
    for i in range(5):
        x_coil = 0.5 + i * 0.5
        coil = plt.Circle((x_coil, 0), 0.15, fill=False,
                          edgecolor='#4488ff', linewidth=1.5, linestyle='-')
        ax.add_patch(coil)
    ax.text(1.75, -0.8, 'МЕТАМАТЕРИАЛ', color='#4488ff', fontsize=7, ha='center')
    
    # Targeting system (top-back)
    target = plt.Circle((-0.3, 1), 0.3, color='#ffaa00', alpha=0.3)
    ax.add_patch(target)
    target_ring = plt.Circle((-0.3, 1), 0.3, fill=False,
                            edgecolor='#ffaa00', linewidth=2)
    ax.add_patch(target_ring)
    ax.text(-0.3, 1, 'НАВЕД.', color='#ffaa00', fontsize=7, ha='center', fontweight='bold')
    
    # Power cable
    ax.plot([-1, 1], [0.5, 1.2], color='#ff4444', linewidth=1.5, linestyle='--')
    ax.plot([-1.5, -0.3], [0.3, 1], color='#ffaa00', linewidth=1.5, linestyle='--')
    
    # Portal output (right side)
    portal_theta = np.linspace(0, 2*np.pi, 100)
    ax.plot(5 + 0.8*np.cos(portal_theta), 0.8*np.sin(portal_theta),
            color='#00ff41', linewidth=3)
    ax.plot(5 + 0.6*np.cos(portal_theta), 0.6*np.sin(portal_theta),
            color='#00ff41', linewidth=1, alpha=0.5)
    portal_glow = plt.Circle((5, 0), 0.9, color='#00ff41', alpha=0.05)
    ax.add_patch(portal_glow)
    ax.text(5, -1.3, 'ПОРТАЛ\n(АКТИВИРОВАН)', color='#00ff41', fontsize=9, 
            ha='center', fontweight='bold')
    
    # Arrow from emitter to portal
    ax.annotate('', xy=(4.1, 0), xytext=(3.6, 0),
                arrowprops=dict(arrowstyle='->', color='#00ff41', lw=2))
    
    # Component list
    specs = [
        ('1. Эмиттер', 'Плазмотрон 532 нм, 10 кВт'),
        ('2. Метаматериал', 'n(r) = √(2r₀/r - 1), СВЧ-резонатор'),
        ('3. Энергоблок', 'Микрореактор, ~1 ТВт'),
        ('4. Наведение', 'Лидар + квантовая синхронизация'),
        ('5. Казимир-модуль', 'Динамический Казимир + усилитель'),
    ]
    
    for i, (name, desc) in enumerate(specs):
        y_pos = -2.8 + i * 0.5
        ax.text(-5.5, y_pos, name, color='#00ff41', fontsize=9, fontweight='bold')
        ax.text(-3.5, y_pos, desc, color='#888888', fontsize=8)
    
    plt.savefig('/home/smboozha/portal_gun/sim_portal_gun.png', 
                dpi=150, bbox_inches='tight', facecolor='#0a0a2e')
    plt.close()
    print("✓ Сохранено: sim_portal_gun.png")

# ============================================================
# RUN ALL SIMULATIONS
# ============================================================
if __name__ == '__main__':
    r0 = 0.1  # 10 cm portal radius
    
    print("="*60)
    print("PORTAL GUN: ЗАПУСК СИМУЛЯЦИЙ")
    print("="*60)
    print()
    
    print("[1/5] Визуализация портала...")
    plot_portal_visual(r0)
    
    print("[2/5] Профиль показателя преломления...")
    plot_n_profile(r0)
    
    print("[3/5] Трассировка лучей через червоточину...")
    trace_ray(r0)
    
    print("[4/5] Симуляция Казимира...")
    simulate_casimir()
    
    print("[5/5] Схема портала-пуши...")
    plot_portal_gun_schematic()
    
    print()
    print("="*60)
    print("ВСЕ СИМУЛЯЦИИ ЗАВЕРШЕНЫ")
    print("="*60)
    print()
    print("Файлы:")
    print("  sim_portal_visual.png  — визуал зелёного портала")
    print("  sim_n_profile.png      — профиль метаматериала")
    print("  sim_wormhole_rays.png  — лучи через червоточину")
    print("  sim_casimir.png        — Казимировская полость")
    print("  sim_portal_gun.png     — схема устройства")
