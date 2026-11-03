"""All main figures for the source-attributed exposure model (3 November 2026), in the 20_publication_maps.py house style.

 1 study area and sources       2 model workflow              3 exposure under two products
 4 fire-attributable exposure   5 industrial and urban excess 6 seasonal decomposition
 7 district source profiles     8 product dependence          9 robustness of the fire effect
10 Manhize comparison          11 control priorities
Usage: python analysis/35_model_figures.py   (after 34_source_model.py)
"""
from pathlib import Path
src = (Path(__file__).parent / '20_publication_maps.py').read_text(encoding='utf-8')
exec(compile(src[:src.index('# Locator and population on exposure support.')], '20_publication_maps_helpers', 'exec'))
import geopandas as gpd
from matplotlib.patches import FancyBboxPatch, Circle
from matplotlib.lines import Line2D
from scipy.stats import spearmanr

SRC = ROOT / 'data' / 'sources'
M = np.load(RES / 'model_cells.npz')
dy = pd.read_csv(RES / 'model_district_year.csv'); summ = pd.read_csv(RES / 'model_district_summary.csv', index_col=0)
coef = pd.read_csv(RES / 'model_fire_coefficients.csv'); mono = pd.read_csv(RES / 'model_fire_monthly_products.csv')
indx = pd.read_csv(RES / 'model_industrial_excess.csv'); daily_p = pd.read_csv(RES / 'model_daily_province.csv', parse_dates=['date'])
J = json.loads((RES / 'source_signatures.json').read_text()); sites = J['hotspots']
its = pd.read_csv(RES / 'source_its.csv')
with rasterio.open(PUB / 'ghap_annual_2022.tif') as s_: GT = s_.transform; GS = s_.shape
GEXT = (GT.c, GT.c + GT.a * GS[1], GT.f + GT.e * GS[0], GT.f)
DIST = list(bounds)
FIRE, IND, BG = '#d9822b', '#7b3294', '#9fb3bf'
SITE_STYLE = {'Kwekwe': '^', 'Redcliff': '^', 'Gweru': 'o', 'Unki': 's', 'Mimosa': 's', 'Manhize': '^'}

def show(ax, arr, cmap, norm, ext=GEXT, clip=True):
    im = ax.imshow(arr, extent=ext, origin='upper', interpolation='nearest', cmap=cmap, norm=norm, zorder=1)
    if clip: im.set_clip_path(PathPatch(region_path(region), transform=ax.transData))
    return im

def mark_sites(ax, names=('Kwekwe', 'Redcliff', 'Gweru', 'Unki', 'Mimosa'), labels=True, ms=6):
    halo = [pe.withStroke(linewidth=2.4, foreground='white', alpha=.95)]
    off = {'Kwekwe': (.07, .05), 'Redcliff': (-.42, -.08), 'Gweru': (-.36, .05), 'Unki': (.07, -.04), 'Mimosa': (.07, -.04), 'Manhize': (-.2, .1)}
    for n in names:
        lo, la = sites[n]
        ax.plot(lo, la, marker=SITE_STYLE[n], ms=ms, mfc='white' if n == 'Gweru' else '#2b2b2b', mec='#2b2b2b', mew=1, zorder=9)
        if labels: ax.text(lo + off[n][0], la + off[n][1], n, fontsize=9, color='#172c38', zorder=9, path_effects=halo)

def bar_style(ax, letter, title=None):
    ax.set_title(letter, loc='left', fontweight='bold', fontsize=12, pad=10)
    if title: ax.text(.06, 1.035, title, transform=ax.transAxes, fontsize=10, va='bottom')
    ax.grid(axis='x', color='#e6eaec', lw=.6); ax.set_axisbelow(True)

# ---------------------------------------------------------------- Figure 1: study area and sources
fig = plt.figure(figsize=(WIDTH, 120 / 25.4))
loc = fig.add_axes([.02, .45, .22, .45])
zw = gpd.read_file(PUB / 'zimbabwe_adm1.geojson')
af = gpd.read_file(PUB / 'natural_earth_110m_countries.geojson')
af.plot(ax=loc, facecolor='#eef1f2', edgecolor='white', lw=.4)
zw.plot(ax=loc, facecolor='#d5dde1', edgecolor='white', lw=.5)
gpd.GeoSeries([region], crs=4326).plot(ax=loc, facecolor=navy, edgecolor=navy)
loc.set_xlim(24.5, 34.5); loc.set_ylim(-23.5, -14.5); loc.set_aspect(1 / np.cos(np.radians(19))); loc.set_xticks([]); loc.set_yticks([])
for sp in loc.spines.values(): sp.set_visible(True); sp.set_color('#9DA9AF')
loc.set_title('a', loc='left', fontweight='bold', fontsize=12, pad=10)
loc.text(29.8, -21.4, 'Zimbabwe', fontsize=9, color='#35424B', ha='center')
ax = panelmap(fig, [.26, .16, .35, .76], 'b')
popd = M['pop2022'] / (1.11 * 1.05)
im = show(ax, np.where(popd > .5, popd, np.nan), 'magma_r', LogNorm(1, 3000)); district_numbers(ax); mark_sites(ax, labels=False)
key(fig, im, [.28, .1, .31, .025], 'Residents per km² (2022)')
ax = panelmap(fig, [.635, .16, .35, .76], 'c')
fs = None
for y in range(2017, 2023):
    for m in range(1, 13):
        with rasterio.open(SRC / f'fire_grid_{y}_{m:02d}.tif') as s_: a = s_.read().sum(0)
        fs = a if fs is None else fs + a
fd = fs * 28.6 / 6 / (5.27 * 5.53) * 100
im = show(ax, np.where(fd > .5, fd, np.nan), 'YlOrBr', LogNorm(1, 200), ext=(27.0, 32.5, -21.9, -16.4))
mark_sites(ax)
key(fig, im, [.655, .1, .31, .025], 'Fire detections per 100 km² per year')
nax = fig.add_axes([.94, .80, .04, .1]); north(nax)
scale_bar(fig, ax, .64, .22, length=100)
hand = [Line2D([], [], marker='^', ls='', mfc='#2b2b2b', mec='#2b2b2b', label='Smelting and steel'),
        Line2D([], [], marker='s', ls='', mfc='#2b2b2b', mec='#2b2b2b', label='Platinum mine'),
        Line2D([], [], marker='o', ls='', mfc='white', mec='#2b2b2b', label='City')]
fig.legend(handles=hand, loc='upper left', bbox_to_anchor=(.02, .40), frameon=False, fontsize=9)
add_key_band(fig); publish(fig, 'Model_Figure_1_study_area')

# ---------------------------------------------------------------- Figure 2: model workflow
fig = plt.figure(figsize=(WIDTH, 95 / 25.4)); ax = fig.add_axes([0, 0, 1, 1]); ax.axis('off'); ax.set_xlim(0, 100); ax.set_ylim(2, 57)
def box(x, y, w, h, title, body, fc):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.4,rounding_size=1.2', fc=fc, ec='#9DA9AF', lw=.8))
    ax.text(x + w / 2, y + h - 2.2, title, ha='center', va='top', fontsize=10, fontweight='bold', color='#172c38')
    ax.text(x + w / 2, y + h - 6.6, body, ha='center', va='top', fontsize=9, color='#35424B', linespacing=1.35)
def arrow(x0, y0, x1, y1):
    ax.annotate('', xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle='-|>', color=navy, lw=1.2))
box(1, 34, 22, 21, 'Exposure', 'GHAP daily PM$_{2.5}$\nACAG annual and\nmonthly PM$_{2.5}$\nWorldPop residents', '#eef4f7')
box(1, 5, 22, 21, 'Sources', 'MODIS active fires\nERA5-Land rain, wind\nTROPOMI NO$_2$\nOpenStreetMap sites', '#fbf1e6')
box(28, 34, 22, 21, 'Fire model', 'district-day panel\nwithin year-month\nln C = β ln(1+F) + ...\nplacebo and checks', '#fbf1e6')
box(28, 5, 22, 21, 'Source detection', 'NO$_2$ excess over\n20 to 50 km ring\nPM$_{2.5}$ local excess\n(10 km), two products', '#f3ecf6')
box(55, 19.5, 20, 21, 'Decomposition', 'every cell and day\nC = fire\n+ industrial or urban\n+ background', '#eef4f7')
box(80, 19.5, 19, 21, 'Attribution', 'exposure above\n15 µg m$^{-3}$ shared\nby source (Shapley)\nper district', '#eaf3ee')
arrow(23.5, 44.5, 27.5, 44.5); arrow(23.5, 15.5, 27.5, 15.5); arrow(23.5, 38, 27.5, 22)
arrow(50.5, 44.5, 54.5, 36); arrow(50.5, 15.5, 54.5, 24); arrow(75.5, 30, 79.5, 30)
publish(fig, 'Model_Figure_2_workflow')

# ---------------------------------------------------------------- Figure 3: exposure under two products
fig = plt.figure(figsize=(WIDTH, 112 / 25.4))
norm = Normalize(9, 15)
for k, (arr, lab) in enumerate(((M['pm_mean'], 'a'), (M['acag_mean'], 'b'))):
    ax = panelmap(fig, [.02 + k * .325, .2, .31, .72], lab); im = show(ax, arr, 'viridis', norm); district_numbers(ax)
    ax.text(.02, .02, 'GHAP' if k == 0 else 'ACAG', transform=ax.transAxes, fontsize=10, fontweight='bold', color='#172c38')
key(fig, im, [.12, .12, .42, .025], 'Mean annual PM$_{2.5}$, 2017 to 2022 (µg m$^{-3}$)')
ax = panelmap(fig, [.67, .2, .31, .72], 'c'); im = show(ax, M['acag_mean'] - M['pm_mean'], 'RdBu_r', TwoSlopeNorm(0, -3, 3)); district_numbers(ax)
key(fig, im, [.69, .12, .27, .025], 'ACAG minus GHAP (µg m$^{-3}$)')
scale_bar(fig, ax, .68, .24, length=100); nax = fig.add_axes([.94, .82, .04, .1]); north(nax)
add_key_band(fig); publish(fig, 'Model_Figure_3_products')

# ---------------------------------------------------------------- Figure 4: fire-attributable exposure
fig = plt.figure(figsize=(WIDTH, 120 / 25.4))
ax = panelmap(fig, [.04, .2, .44, .72], 'a'); im = show(ax, M['fire_inc_son'], 'YlOrBr', Normalize(0, 5)); district_numbers(ax)
key(fig, im, [.08, .12, .36, .025], 'Fire-attributable PM$_{2.5}$, Sep to Nov mean (µg m$^{-3}$)')
ax = panelmap(fig, [.52, .2, .44, .72], 'b'); im = show(ax, 100 * M['fire_share'], 'OrRd', Normalize(0, 25)); district_numbers(ax)
key(fig, im, [.56, .12, .36, .025], 'Share of days above 15 µg m$^{-3}$ caused by fires (%)')
scale_bar(fig, ax, .55, .26, length=100); nax = fig.add_axes([.94, .82, .04, .1]); north(nax)
add_key_band(fig); publish(fig, 'Model_Figure_4_fire')

# ---------------------------------------------------------------- Figure 5: industrial and urban local excess
fig = plt.figure(figsize=(WIDTH, 100 / 25.4))
ax = fig.add_axes([.04, .22, .46, .72]); ax.set_facecolor('#f4f6f6')
win = (29.45, 30.25, -19.75, -18.6)
im = ax.imshow(np.where(M['ind_mean'] > .05, M['ind_mean'], np.nan), extent=GEXT, origin='upper', interpolation='nearest',
               cmap='Purples', norm=Normalize(0, 4), zorder=1)
bg = ax.imshow(M['pm_mean'], extent=GEXT, origin='upper', interpolation='nearest', cmap='Greys', norm=Normalize(8, 18), alpha=.25, zorder=0)
for n, g in bounds.items(): outline(ax, g, lw=.6)
for n in ('Kwekwe', 'Redcliff', 'Gweru'):
    lo, la = sites[n]; ax.add_patch(Circle((lo, la), 10 / 111.32 / np.cos(np.radians(19)), fill=False, ec='#2b2b2b', lw=.8, ls=(0, (3, 2)), zorder=6))
mark_sites(ax, names=('Kwekwe', 'Redcliff', 'Gweru', 'Unki'), labels=False)
halo = [pe.withStroke(linewidth=2.4, foreground='white', alpha=.95)]
for n, (dx, dy_) in {'Kwekwe': (.1, .03), 'Redcliff': (-.27, -.04), 'Gweru': (.1, -.02), 'Unki': (.04, .04)}.items():
    ax.text(sites[n][0] + dx, sites[n][1] + dy_, n, fontsize=9, color='#172c38', zorder=9, path_effects=halo)
ax.set_xlim(win[:2]); ax.set_ylim(win[2:]); ax.set_aspect(1 / np.cos(np.radians(19))); ax.set_xticks([]); ax.set_yticks([])
ax.set_title('a', loc='left', fontweight='bold', fontsize=12, pad=10)
key(fig, im, [.1, .13, .34, .025], 'Local PM$_{2.5}$ excess, GHAP 2017 to 2022 (µg m$^{-3}$)')
ix = indx.groupby('site')[['ghap_excess', 'acag_excess']].agg(['mean', 'min', 'max'])
order = ['Gweru', 'Kwekwe', 'Redcliff', 'Mimosa', 'Unki']
ax = fig.add_axes([.62, .18, .35, .7]); yv = np.arange(len(order))[::-1]
for off_, (prod, c) in zip((.15, -.15), (('ghap_excess', navy), ('acag_excess', teal))):
    r = ix.loc[order, prod]
    ax.hlines(yv + off_, r['min'], r['max'], color=c, lw=1.4); ax.plot(r['mean'], yv + off_, 'o', color=c, ms=5, label=prod.split('_')[0].upper())
ax.axvline(0, color='#9DA9AF', lw=.8); ax.set_yticks(yv); ax.set_yticklabels(order)
ax.set_xlabel('Local excess over 20 to 50 km ring\n(µg m$^{-3}$; mean and range of years)'); ax.legend(frameon=False, loc='lower right')
bar_style(ax, 'b')
publish(fig, 'Model_Figure_5_industry')

# ---------------------------------------------------------------- Figure 6: seasonal decomposition
fig = plt.figure(figsize=(WIDTH, 120 / 25.4))
dp = daily_p.set_index('date').rolling(7, center=True, min_periods=4).mean()
ax = fig.add_axes([.08, .56, .9, .36])
ax.stackplot(dp.index, dp.bg, dp.ind, dp.fire, colors=[BG, IND, FIRE], lw=0, labels=['Background', 'Industrial and urban', 'Fires within 50 km'])
ax.axhline(15, color='#8c2d04', lw=.8, ls=(0, (4, 2))); ax.text(dp.index[5], 15.5, 'WHO 24 h guideline', fontsize=9, color='#8c2d04')
ax.set_ylabel('PM$_{2.5}$ (µg m$^{-3}$), 7-day mean'); ax.set_xlim(dp.index[0], dp.index[-1])
ax.legend(frameon=False, ncol=3, loc='upper left', fontsize=9); ax.set_title('a', loc='left', fontweight='bold', fontsize=12, pad=10)
mc = daily_p.assign(m=daily_p.date.dt.month).groupby('m')[['bg', 'ind', 'fire']].mean()
ax = fig.add_axes([.08, .1, .55, .34]); b0 = np.zeros(12)
for col, c in (('bg', BG), ('ind', IND), ('fire', FIRE)):
    ax.bar(mc.index, mc[col], bottom=b0, color=c, width=.75); b0 += mc[col].to_numpy()
ax.set_xticks(range(1, 13)); ax.set_xticklabels(list('JFMAMJJASOND')); ax.set_ylabel('PM$_{2.5}$ (µg m$^{-3}$)')
ax.set_title('b', loc='left', fontweight='bold', fontsize=12, pad=10)
sh = 100 * mc.fire / mc.sum(1)
ax = fig.add_axes([.72, .1, .26, .34]); ax.plot(sh.index, sh, '-o', color=FIRE, ms=4)
ax.set_xticks(range(1, 13)); ax.set_xticklabels(list('JFMAMJJASOND')); ax.set_ylabel('Fire share of PM$_{2.5}$ (%)')
ax.set_title('c', loc='left', fontweight='bold', fontsize=12, pad=10)
publish(fig, 'Model_Figure_6_seasonal')

# ---------------------------------------------------------------- Figure 7: district source profiles
fig = plt.figure(figsize=(WIDTH, 100 / 25.4))
s7 = summ.drop('Midlands Province').sort_values('pd_total')
ax = fig.add_axes([.14, .16, .38, .72]); yv = np.arange(len(s7)); left = np.zeros(len(s7))
for col, c, lab in (('pd_other', BG, 'Background'), ('pd_industry', IND, 'Industrial and urban'), ('pd_fire', FIRE, 'Fires')):
    ax.barh(yv, s7[col] / 1e6, left=left, color=c, height=.7, label=lab); left += (s7[col] / 1e6).to_numpy()
for yi, (n, r) in zip(yv, s7.iterrows()):
    ax.text(left[yi] + 3, yi, f"{r.fire_share + r.ind_share:.0f}%", va='center', fontsize=9, color='#35424B')
ax.set_yticks(yv); ax.set_yticklabels(s7.index); ax.set_xlabel('Population-days above 15 µg m$^{-3}$, 2017 to 2022 (millions)')
ax.legend(frameon=False, loc='lower right', fontsize=9); bar_style(ax, 'a'); ax.set_xlim(0, 260)
s8 = summ.drop('Midlands Province').sort_values('pm')
ax = fig.add_axes([.64, .16, .34, .72]); left = np.zeros(len(s8))
for col, c in (('background', BG), ('industry', IND), ('fire', FIRE)):
    ax.barh(np.arange(len(s8)), s8[col], left=left, color=c, height=.7); left += s8[col].to_numpy()
ax.set_yticks(np.arange(len(s8))); ax.set_yticklabels(s8.index); ax.set_xlabel('Mean PM$_{2.5}$, 2017 to 2022 (µg m$^{-3}$)')
ax.axvline(5, color='#8c2d04', lw=.8, ls=(0, (4, 2))); bar_style(ax, 'b')
publish(fig, 'Model_Figure_7_districts')

# ---------------------------------------------------------------- Figure 8: product dependence
fig = plt.figure(figsize=(WIDTH, 92 / 25.4))
ax = fig.add_axes([.15, .2, .2, .66])
rows8 = [('GHAP daily', coef.iloc[0].pct_per_doubling, coef.iloc[0].lo, coef.iloc[0].hi)] + \
        [(f"{r['product']} monthly", r.pct_per_doubling, r.lo, r.hi) for _, r in mono.iterrows()]
for yi, (lab, e, lo, hi) in enumerate(rows8[::-1]):
    ax.plot([lo, hi], [yi, yi], color=navy, lw=1.4); ax.plot(e, yi, 'o', color=FIRE if 'daily' in lab else navy, ms=5)
ax.axvline(0, color='#9DA9AF', lw=.8); ax.set_yticks(range(len(rows8))); ax.set_yticklabels([r[0] for r in rows8[::-1]])
ax.set_xlabel('PM$_{2.5}$ change per doubling of fires (%)'); bar_style(ax, 'a')
dist_ = summ.drop('Midlands Province')
z22 = np.load(RES / 'population_weights_primary_2022.npz')
acag_pw = {}
for d in dist_.index:
    w = z22[d.replace(' ', '_')].reshape(GS); a = M['acag_mean']; ok = np.isfinite(a) & (w > 0)
    acag_pw[d] = (w[ok] * a[ok]).sum() / w[ok].sum()
dist_ = dist_.assign(acag=pd.Series(acag_pw))
dist_['acag_ind'] = [indx[indx.site.isin(['Kwekwe', 'Redcliff'])].acag_excess.mean() * dist_.loc[d, 'industry'] /
                     max(indx[indx.site.isin(['Kwekwe', 'Redcliff'])].ghap_excess.mean(), 1e-9) if d == 'Kwekwe' else
                     (indx[indx.site == 'Gweru'].acag_excess.mean() / indx[indx.site == 'Gweru'].ghap_excess.mean() * dist_.loc[d, 'industry']
                      if d == 'Gweru' else 0) for d in dist_.index]
dist_['acag_bg'] = dist_.acag - dist_.fire * dist_.acag / dist_.pm - dist_.acag_ind
ranks = pd.DataFrame({'GHAP total': dist_.pm.rank(ascending=False), 'ACAG total': dist_.acag.rank(ascending=False),
                      'GHAP background': dist_.background.rank(ascending=False), 'ACAG background': dist_.acag_bg.rank(ascending=False)})
ax = fig.add_axes([.43, .2, .2, .66])
for d, r in ranks.iterrows():
    hl = d in ('Gweru', 'Kwekwe'); c = IND if hl else '#9DA9AF'
    ax.plot([0, 1], [r['GHAP total'], r['ACAG total']], '-o', color=c, lw=1.6 if hl else 1, ms=4)
    ax.text(-.06, r['GHAP total'], d, ha='right', va='center', fontsize=9, color='#172c38' if hl else '#54636d')
ax.set_xticks([0, 1]); ax.set_xticklabels(['GHAP', 'ACAG']); ax.invert_yaxis(); ax.set_yticks([])
ax.set_xlim(-.05, 1.1); ax.set_title('b', loc='left', fontweight='bold', fontsize=12, pad=10)
ax.text(.6, 1.035, f"Total, ρ ={spearmanr(ranks['GHAP total'], ranks['ACAG total'])[0]:.2f}", transform=ax.transAxes, ha='center', fontsize=10)
ax = fig.add_axes([.78, .2, .2, .66])
for d, r in ranks.iterrows():
    hl = d in ('Gweru', 'Kwekwe'); c = IND if hl else '#9DA9AF'
    ax.plot([0, 1], [r['GHAP background'], r['ACAG background']], '-o', color=c, lw=1.6 if hl else 1, ms=4)
    ax.text(-.06, r['GHAP background'], d, ha='right', va='center', fontsize=9, color='#172c38' if hl else '#54636d')
ax.set_xticks([0, 1]); ax.set_xticklabels(['GHAP', 'ACAG']); ax.invert_yaxis(); ax.set_yticks([])
ax.set_xlim(-.05, 1.1); ax.set_title('c', loc='left', fontweight='bold', fontsize=12, pad=10)
ax.text(.62, 1.035, f"Background, ρ ={spearmanr(ranks['GHAP background'], ranks['ACAG background'])[0]:.2f}", transform=ax.transAxes, ha='center', fontsize=10)
ranks.to_csv(RES / 'model_district_ranks.csv')
publish(fig, 'Model_Figure_8_products')

# ---------------------------------------------------------------- Figure 9: robustness of the fire effect
fig = plt.figure(figsize=(WIDTH, 120 / 25.4)); ax = fig.add_axes([.32, .1, .64, .84])
c9 = coef[coef.model != 'Main, with later fires held fixed'].reset_index(drop=True)
yv = np.arange(len(c9))[::-1]
for yi, (_, r) in zip(yv, c9.iterrows()):
    col = '#8c2d04' if r.model.startswith('Placebo') else (FIRE if r.model.startswith('Main') else navy)
    ax.plot([r.lo, r.hi], [yi, yi], color=col, lw=1.4); ax.plot(r.pct_per_doubling, yi, 'o', color=col, ms=5)
ax.axvline(0, color='#9DA9AF', lw=.8); ax.axvline(c9.iloc[0].pct_per_doubling, color=FIRE, lw=.7, ls=(0, (3, 2)))
ax.set_yticks(yv); ax.set_yticklabels(c9.model); ax.set_xlabel('PM$_{2.5}$ change per doubling of fires (%), 95% CI')
ax.grid(axis='x', color='#e6eaec', lw=.6); ax.set_axisbelow(True)
publish(fig, 'Model_Figure_9_robustness')

# ---------------------------------------------------------------- Figure 10: Manhize comparison
ms = pd.read_csv(SRC / 'monthly_sites.csv'); ms['date'] = pd.to_datetime(ms.month + '-01'); ms['ex'] = ms.no2_near - ms.no2_ring
fig = plt.figure(figsize=(WIDTH, 85 / 25.4))
ax = fig.add_axes([.08, .17, .52, .7])
cmp_ = ms[ms.site != 'Manhize'].groupby('date').ex.mean().rolling(3, center=True, min_periods=2).mean()
mz = ms[ms.site == 'Manhize'].set_index('date').ex.rolling(3, center=True, min_periods=2).mean()
ax.plot(cmp_.index, cmp_, color=navy, lw=1.4, label='Five Midlands sites (mean)'); ax.plot(mz.index, mz, color=FIRE, lw=1.4, label='Manhize')
ax.axvline(pd.Timestamp('2024-06-12'), color='#2b2b2b', lw=.8, ls=(0, (4, 2)))
ax.text(pd.Timestamp('2024-04-01'), -.85, 'Blast furnace\nfired, June 2024', fontsize=9, va='bottom', ha='right')
ax.set_ylabel('Local NO$_2$ excess (µmol m$^{-2}$), 3-month mean'); ax.legend(frameon=False, loc='upper left', fontsize=9)
ax.set_title('a', loc='left', fontweight='bold', fontsize=12, pad=10)
g = its[its.gas == 'no2'].set_index('site').loc[['Manhize', 'Gweru', 'Kwekwe', 'Unki', 'Redcliff', 'Mimosa']]
ax = fig.add_axes([.72, .17, .26, .7]); yv = np.arange(len(g))[::-1]
for yi, (n, r) in zip(yv, g.iterrows()):
    c = FIRE if n == 'Manhize' else navy; ax.plot([r.ci_low, r.ci_high], [yi, yi], color=c, lw=1.4); ax.plot(r.effect, yi, 'o', color=c, ms=5)
ax.axvline(0, color='#9DA9AF', lw=.8); ax.set_yticks(yv); ax.set_yticklabels(g.index); ax.set_xlabel('Change after June 2024')
bar_style(ax, 'b')
publish(fig, 'Model_Figure_10_manhize')

# ---------------------------------------------------------------- Figure 11: control priorities
fig = plt.figure(figsize=(WIDTH, 112 / 25.4))
ax = panelmap(fig, [.03, .1, .5, .82], 'a')
def cls(r):
    f_, i_ = r.fire_share >= 12, r.ind_share >= 5
    return 'Fires and local sources' if f_ and i_ else 'Fires' if f_ else 'Local sources' if i_ else 'Regional background'
CC = {'Fires and local sources': '#b35806', 'Fires': '#f1a340', 'Local sources': IND, 'Regional background': '#d8dee2'}
s11 = summ.drop('Midlands Province').copy(); s11['class'] = s11.apply(cls, axis=1)
for n, g_ in bounds.items():
    outline(ax, g_, lw=.6, fill=CC[s11.loc[n, 'class']])
district_numbers(ax); mark_sites(ax, names=('Kwekwe', 'Redcliff', 'Gweru'), labels=False)
fig.legend(handles=[Patch(fc=c, ec='#9DA9AF', label=k) for k, c in CC.items() if k in set(s11['class'])], loc='lower left',
           bbox_to_anchor=(.05, .0), ncol=3, frameon=False, fontsize=9)
ax = fig.add_axes([.64, .2, .34, .7]); yv = np.arange(len(s11))
s11 = s11.sort_values('fire_share')
ax.barh(yv + .18, s11.fire_share, height=.34, color=FIRE, label='Fires')
ax.barh(yv - .18, s11.ind_share, height=.34, color=IND, label='Industrial and urban')
ax.set_yticks(yv); ax.set_yticklabels(s11.index); ax.set_xlabel('Share of exceedance population-days (%)')
ax.legend(frameon=False, loc='lower right', fontsize=9); bar_style(ax, 'b')
s11.to_csv(RES / 'model_priorities.csv')
add_key_band(fig); publish(fig, 'Model_Figure_11_priorities')
