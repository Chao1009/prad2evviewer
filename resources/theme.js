// theme.js — runtime colour theme for the web monitor.
//
// All CSS colours are defined as custom properties in viewer.css: classic in
// the default :root scope, dark/light under :root[data-theme=...]. This module
// reads those properties at load time and on every switch, exposes them as
// the THEME object to the rest of the JS, and provides Plotly-layout helpers
// so plots pick up the active palette.
//
// Public API:
//   THEME            — object with bg/canvas/text/.../accent/danger/... keys
//   currentTheme()   — returns 'dark', 'light' or 'classic'
//   setTheme(name)   — flips data-theme, persists to localStorage, notifies
//   toggleTheme()    — cycle to the next theme
//   onThemeChange(fn) — register a callback (fn(newName))
//   plotlyLayout()   — Plotly layout skeleton for the active theme
//   plotlyRelayout(divId) — apply the active theme layout to an existing plot

'use strict';

const THEME_STORAGE_KEY = 'prad2.theme';
const THEME_NAMES = ['dark', 'light', 'classic'];
const THEME = {};
const _themeListeners = [];

// --theme-* custom properties mirrored into THEME; each key is the camelCase
// form of the token name (select-border -> THEME.selectBorder).
const THEME_TOKENS = [
    'bg', 'canvas', 'alt-base',
    'border', 'grid', 'text', 'text-dim', 'text-muted',
    'accent',
    'success', 'warn', 'danger', 'highlight', 'no-data', 'select-border',
    'overlay', 'cut-shade',
];

function refreshTheme(){
    const cs = getComputedStyle(document.documentElement);
    for(const t of THEME_TOKENS)
        THEME[t.replace(/-(\w)/g, (_, c) => c.toUpperCase())] =
            cs.getPropertyValue('--theme-' + t).trim();
}

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

// Plotly helpers

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
