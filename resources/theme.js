// theme.js — runtime Apple-inspired theme for the web monitor.
//
// All CSS colours are defined as custom properties in viewer.css under
// :root[data-theme="dark"] / [data-theme="light"]. This module reads those
// properties at load time and on every switch, exposes them as the THEME
// object to the rest of the JS, and provides Plotly-layout helpers so plots
// pick up the active palette.
//
// Public API:
//   THEME            — object with bg/canvas/panel/.../accent/danger/... keys
//   currentTheme()   — returns 'dark' or 'light'
//   setTheme(name)   — flips data-theme, persists to localStorage, notifies
//   onThemeChange(fn) — register a callback (fn(newName))
//   plotlyLayout()   — Plotly layout skeleton for the active theme
//   plotlyRelayout(divId) — apply the active theme layout to an existing plot

'use strict';

const THEME_STORAGE_KEY = 'prad2.theme';
const THEME_NAMES = ['dark', 'light', 'classic'];
const THEME = {};
const _themeListeners = [];

// Fallbacks matching the "classic" palette in viewer.css (the default) —
// used only if the CSS hasn't finished loading yet when this script runs.
// Once the stylesheet applies, refreshTheme() picks up the real values.
const _FALLBACK = {
    '--theme-bg':            '#1a1a2e',
    '--theme-canvas':        '#11112a',
    '--theme-panel':         '#16213e',
    '--theme-button':        '#0f3460',
    '--theme-button-hover':  '#1f2d4f',
    '--theme-alt-base':      '#0a0a18',
    '--theme-tooltip':       'rgba(0,0,0,0.92)',
    '--theme-border':        '#0f3460',
    '--theme-grid':          '#333333',
    '--theme-text':          '#e0e0e0',
    '--theme-text-strong':   '#ffffff',
    '--theme-text-dim':      '#888888',
    '--theme-text-muted':    '#555555',
    '--theme-accent':        '#00b4d8',
    '--theme-accent-strong': '#0074d9',
    '--theme-accent-border': '#00b4d8',
    '--theme-on-accent':     '#000000',
    '--theme-success':       '#51cf66',
    '--theme-warn':          '#ffa500',
    '--theme-danger':        '#ff6b6b',
    '--theme-highlight':     '#ff922b',
    '--theme-no-data':       '#1a1a2e',
    '--theme-select-border': '#ffffff',
    '--theme-overlay':       'rgba(0,0,0,0.7)',
    '--theme-overlay-light': 'rgba(0,0,0,0.35)',
    '--theme-shadow':        'rgba(0,0,0,0.5)',
};

function _readCssVar(name){
    const v = getComputedStyle(document.documentElement)
        .getPropertyValue(name).trim();
    return v || _FALLBACK[name] || '';
}

function refreshTheme(){
    const keys = [
        ['bg',           '--theme-bg'],
        ['canvas',       '--theme-canvas'],
        ['panel',        '--theme-panel'],
        ['button',       '--theme-button'],
        ['buttonHover',  '--theme-button-hover'],
        ['altBase',      '--theme-alt-base'],
        ['tooltip',      '--theme-tooltip'],
        ['border',       '--theme-border'],
        ['grid',         '--theme-grid'],
        ['text',         '--theme-text'],
        ['textStrong',   '--theme-text-strong'],
        ['textDim',      '--theme-text-dim'],
        ['textMuted',    '--theme-text-muted'],
        ['accent',       '--theme-accent'],
        ['accentStrong', '--theme-accent-strong'],
        ['accentBorder', '--theme-accent-border'],
        ['onAccent',     '--theme-on-accent'],
        ['success',      '--theme-success'],
        ['warn',         '--theme-warn'],
        ['danger',       '--theme-danger'],
        ['highlight',    '--theme-highlight'],
        ['noData',       '--theme-no-data'],
        ['selectBorder', '--theme-select-border'],
        ['overlay',      '--theme-overlay'],
        ['overlayLight', '--theme-overlay-light'],
        ['shadow',       '--theme-shadow'],
    ];
    for(const [k, css] of keys) THEME[k] = _readCssVar(css);
}

function availableThemes(){ return THEME_NAMES.slice(); }

function currentTheme(){
    return document.documentElement.dataset.theme || 'classic';
}

function setTheme(name){
    if(!THEME_NAMES.includes(name)) return;
    if(currentTheme() === name) return;
    document.documentElement.dataset.theme = name;
    try { localStorage.setItem(THEME_STORAGE_KEY, name); } catch(e){}
    refreshTheme();
    for(const fn of _themeListeners) { try { fn(name); } catch(e){ console.error(e); } }
}

// Cycle dark → light → classic → dark …
function toggleTheme(){
    const i = THEME_NAMES.indexOf(currentTheme());
    const next = THEME_NAMES[(i + 1 + THEME_NAMES.length) % THEME_NAMES.length];
    setTheme(next);
}

function onThemeChange(fn){ _themeListeners.push(fn); }

// Initialise before first paint.
(function initTheme(){
    let saved = null;
    try { saved = localStorage.getItem(THEME_STORAGE_KEY); } catch(e){}
    document.documentElement.dataset.theme =
        THEME_NAMES.includes(saved) ? saved : 'classic';
    refreshTheme();
})();

// -------------------------------------------------------------------------
// Plotly helpers
// -------------------------------------------------------------------------

// Base layout skeleton matching the active theme. Callers spread this into
// their layout and then add their own title / margin / axis titles.
function plotlyLayout(){
    return {
        paper_bgcolor: THEME.bg,
        plot_bgcolor:  THEME.canvas,
        font: { family: 'Consolas,monospace', size: 10, color: THEME.textDim },
        margin: { l: 45, r: 10, t: 24, b: 32 },
        xaxis: {
            gridcolor: THEME.grid, zerolinecolor: THEME.border,
            linecolor: THEME.border, tickcolor: THEME.border,
        },
        yaxis: {
            gridcolor: THEME.grid, zerolinecolor: THEME.border,
            linecolor: THEME.border, tickcolor: THEME.border,
        },
    };
}

// Keys to overwrite on an existing plot when the theme flips.
function plotlyThemePatch(){
    return {
        paper_bgcolor: THEME.bg,
        plot_bgcolor:  THEME.canvas,
        'font.color':         THEME.textDim,
        'xaxis.gridcolor':    THEME.grid,
        'xaxis.zerolinecolor':THEME.border,
        'xaxis.linecolor':    THEME.border,
        'xaxis.tickcolor':    THEME.border,
        'yaxis.gridcolor':    THEME.grid,
        'yaxis.zerolinecolor':THEME.border,
        'yaxis.linecolor':    THEME.border,
        'yaxis.tickcolor':    THEME.border,
    };
}

// Re-apply the active theme's chrome to a plot that was drawn earlier.
function plotlyRelayout(divId){
    try {
        if(typeof Plotly !== 'undefined' && document.getElementById(divId)) {
            Plotly.relayout(divId, plotlyThemePatch());
        }
    } catch(e){ /* plot may not be initialised yet */ }
}

// Expose everything on the global window so non-module scripts can use it.
window.THEME            = THEME;
window.availableThemes  = availableThemes;
window.refreshTheme     = refreshTheme;
window.currentTheme     = currentTheme;
window.setTheme         = setTheme;
window.toggleTheme      = toggleTheme;
window.onThemeChange    = onThemeChange;
window.plotlyLayout     = plotlyLayout;
window.plotlyThemePatch = plotlyThemePatch;
window.plotlyRelayout   = plotlyRelayout;
