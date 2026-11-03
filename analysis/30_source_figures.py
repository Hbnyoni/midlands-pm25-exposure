"""Figures 10 and 11 for the GeoHealth resubmission (3 November 2026), in the 20_publication_maps.py house style.

Figure 10  Fires and the September to November exposure peak (a monthly cycle, b within-month association, c high and low
           fire days in September to November).
Figure 11  Industrial and urban source signatures (a TROPOMI NO2 before June 2024 with sites, b change in local NO2 excess
           after June 2024 by site, c GHAP PM2.5 near the Kwekwe industrial belt against its surroundings).
Usage: python analysis/30_source_figures.py   (after 28 and 29)
"""
from pathlib import Path
src = (Path(__file__).parent / '20_publication_maps.py').read_text(encoding='utf-8')
exec(compile(src[:src.index('# Locator and population on exposure support.')], '20_publication_maps_helpers', 'exec'))
import statsmodels.formula.api as smf

SRC = ROOT / 'data' / 'sources'
J = json.loads((RES / 'source_signatures.json').read_text())
its = pd.read_csv(RES / 'source_its.csv')
d = pd.read_csv(RES / 'source_fire_daily.csv', parse_dates=['date'])
FIRE = '#c7772e'

# ---- Figure 10 -------------------------------------------------------------------------------------------------
fig = plt.figure(figsize=(WIDTH, 120 / 25.4))
mo = d.groupby('month').agg(fires=('fires', 'sum'), pm=('pw_pm25', 'mean'))
mo['fires_per_year'] = mo.fires / 6
ax = fig.add_axes([.11, .58, .37, .33])
ax.bar(mo.index, mo.fires_per_year, color=FIRE, alpha=.85, width=.72)
ax.axvspan(8.5, 11.5, color='#f2e6d8', zorder=0)
ax.set_xticks(range(1, 13)); ax.set_xticklabels(list('JFMAMJJASOND')); ax.set_ylabel('Fire pixels per month')
ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v / 1000:g}k' if v else '0'))
ax.set_title('a', loc='left', fontweight='bold', fontsize=12, pad=10)
ax2 = fig.add_axes([.11, .14, .37, .33], sharex=ax)
ax2.plot(mo.index, mo.pm, '-o', color=navy, ms=4, lw=1.6)
ax2.axhline(15, color='#8c2d04', lw=.8, ls=(0, (4, 2))); ax2.text(1.1, 15.3, 'WHO 24 h guideline', fontsize=9, color='#8c2d04')
ax2.axvspan(8.5, 11.5, color='#f2e6d8', zorder=0); ax2.set_ylabel('PM$_{2.5}$ (µg m$^{-3}$)')
ax2.set_xticks(range(1, 13)); ax2.set_xticklabels(list('JFMAMJJASOND'))
dd = d.dropna(subset=['lpm', 'lfire', 'precip_mm', 'wind']).copy()
dd['r_pm'] = smf.ols('lpm ~ precip_mm + wind + C(ym)', dd).fit().resid
dd['r_fire'] = smf.ols('lfire ~ precip_mm + wind + C(ym)', dd).fit().resid
bins = pd.qcut(dd.r_fire, 15, duplicates='drop')
b = dd.groupby(bins, observed=True).agg(x=('r_fire', 'mean'), y=('r_pm', 'mean'), n=('r_pm', 'size'), s=('r_pm', 'sem'))
ax = fig.add_axes([.585, .14, .17, .77])
ax.scatter(dd.r_fire, dd.r_pm, s=3, color=muted, alpha=.12, lw=0)
ax.errorbar(b.x, b.y, yerr=1.96 * b.s, fmt='o', color=navy, ms=4, lw=.8, capsize=0)
xs = np.linspace(dd.r_fire.quantile(.01), dd.r_fire.quantile(.99), 50); ax.plot(xs, J['fires']['b_log_fire'] * xs, color=FIRE, lw=1.6)
ax.set_xlim(xs[0], xs[-1]); ax.set_ylim(-.45, .45)
ax.set_xlabel('Fires, 3-day log anomaly'); ax.set_ylabel('PM$_{2.5}$, log anomaly')
ax.set_title('b', loc='left', fontweight='bold', fontsize=12, pad=10)
f = J['fires']
ax.text(.03, .97, f"+{f['pct_pm_per_doubling']:.1f}% per doubling\n(95% CI {f['pct_ci'][0]:.1f} to {f['pct_ci'][1]:.1f})",
        transform=ax.transAxes, va='top', fontsize=9, color='#202B33')
son = d[d.month.isin([9, 10, 11])].dropna(subset=['pw_pm25', 'fires3'])
hi = son.fires3 > son.fires3.median()
ax = fig.add_axes([.83, .14, .15, .77])
parts_ = ax.violinplot([son[~hi].pw_pm25, son[hi].pw_pm25], positions=[0, 1], widths=.75, showmedians=True, showextrema=False)
for body, c in zip(parts_['bodies'], (muted, FIRE)): body.set_facecolor(c); body.set_alpha(.55); body.set_edgecolor('none')
parts_['cmedians'].set_color(navy)
ax.axhline(15, color='#8c2d04', lw=.8, ls=(0, (4, 2)))
ax.set_xticks([0, 1]); ax.set_xticklabels(['Fewer\nfires', 'More\nfires']); ax.set_ylabel('PM$_{2.5}$ (µg m$^{-3}$), Sep to Nov')
ax.set_title('c', loc='left', fontweight='bold', fontsize=12, pad=10)
publish(fig, 'Figure_10_fires')

# ---- Figure 11 -------------------------------------------------------------------------------------------------
fig = plt.figure(figsize=(WIDTH, 130 / 25.4))
ax = panelmap(fig, [.02, .30, .40, .64], 'a')
with rasterio.open(SRC / 'tropomi_no2_pre.tif') as s_:
    a = s_.read(1).astype(float); bb = s_.bounds
a[a < -9000] = np.nan
im = ax.imshow(a, extent=(bb.left, bb.right, bb.bottom, bb.top), origin='upper', cmap='magma_r', norm=Normalize(12, 19),
               interpolation='nearest', zorder=1)
clip = region.buffer(.3)
im.set_clip_path(PathPatch(region_path(clip), transform=ax.transData))
outline(ax, region, lw=1.1)
cb_ = clip.bounds; ax.set_xlim(cb_[0], cb_[2]); ax.set_ylim(cb_[1], cb_[3])
SITE = {'Manhize': 'Manhize\n(steel, 2024)', 'Kwekwe': 'Kwekwe', 'Redcliff': 'Redcliff', 'Unki': 'Unki', 'Mimosa': 'Mimosa', 'Gweru': 'Gweru'}
OFF = {'Manhize': (.06, -.20), 'Kwekwe': (.08, .05), 'Redcliff': (-.42, -.07), 'Unki': (.08, -.03), 'Mimosa': (.08, -.03), 'Gweru': (-.33, .03)}
halo = [pe.withStroke(linewidth=2.4, foreground='white', alpha=.95)]
for k, (lo, la) in J['hotspots'].items():
    ax.plot(lo, la, marker='^' if k != 'Gweru' else 'o', ms=6.5, mfc='#38a3a5' if k != 'Gweru' else 'white', mec=navy, mew=.9, zorder=8)
    ax.text(lo + OFF[k][0], la + OFF[k][1], SITE[k], fontsize=9, color='#172c38', zorder=8, path_effects=halo)
key(fig, im, [.05, .16, .36, .025], 'Tropospheric NO$_2$, July 2018\nto May 2024 (µmol m$^{-2}$)')
scale_bar(fig, ax, .26, .225, length=100)
nax = fig.add_axes([.37, .78, .04, .1]); north(nax)
g = its[its.gas == 'no2'].set_index('site').loc[['Manhize', 'Gweru', 'Kwekwe', 'Unki', 'Redcliff', 'Mimosa']]
ax = fig.add_axes([.58, .58, .39, .33])
y = np.arange(len(g))[::-1]
for yi, (site, r) in zip(y, g.iterrows()):
    c = FIRE if site == 'Manhize' else navy
    ax.plot([r.ci_low, r.ci_high], [yi, yi], color=c, lw=1.4); ax.plot(r.effect, yi, 'o', color=c, ms=5)
ax.axvline(0, color='#9DA9AF', lw=.8)
ax.set_yticks(y); ax.set_yticklabels(g.index); ax.set_xlabel('Change in NO$_2$ excess after June 2024 (µmol m$^{-2}$)')
ax.set_title('b', loc='left', fontweight='bold', fontsize=12, pad=10)
rows = []
for yr in range(2017, 2023):
    with rasterio.open(PUB / f'ghap_annual_{yr}.tif') as s_:
        aa = s_.read(1).astype(float); t = s_.transform
    aa[(aa < -9000) | (aa <= 0)] = np.nan
    yy, xx = np.mgrid[0:aa.shape[0], 0:aa.shape[1]]; lon = t.c + (xx + .5) * t.a; lat = t.f + (yy + .5) * t.e
    lo0, la0 = J['hotspots']['Kwekwe']
    dist = np.hypot((lon - lo0) * 111.32 * np.cos(np.radians(la0)), (lat - la0) * 110.57)
    rows.append((yr, np.nanmean(aa[dist <= 10]), np.nanmean(aa[(dist >= 20) & (dist <= 50)])))
kk = pd.DataFrame(rows, columns=['year', 'near', 'ring'])
ax = fig.add_axes([.53, .14, .44, .30])
ax.plot(kk.year, kk.near, '-o', color=FIRE, ms=4, lw=1.5, label='Within 10 km of the Kwekwe NO$_2$ maximum')
ax.plot(kk.year, kk.ring, '-o', color=navy, ms=4, lw=1.5, label='20 to 50 km ring')
ax.set_ylabel('GHAP annual PM$_{2.5}$ (µg m$^{-3}$)'); ax.legend(frameon=False, fontsize=9, loc='lower left')
ax.set_title('c', loc='left', fontweight='bold', fontsize=12, pad=10); ax.set_ylim(7, 17)
publish(fig, 'Figure_11_sources')
