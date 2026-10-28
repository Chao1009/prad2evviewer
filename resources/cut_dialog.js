// cut_dialog.js — Waveform-Tab peak filter dialog ("Cut Settings…")
// Called from init() in viewer.js.

function initCutDialog(){
    const $ = id => document.getElementById(id);
    if (!$('cut-backdrop') || !$('cut-dialog')) return;

    // Build quality-bit checkbox lists from histConfig.quality_bits.
    // Always rebuilds — cheap, and avoids subtle bugs when the bit palette
    // changes (e.g. server config reloaded between opens).
    function buildBitList(containerId, set){
        const c = $(containerId);
        if (!c) return;
        const bits = histConfig.quality_bits || [];
        c.innerHTML = '';
        if (!bits.length) {
            const empty = document.createElement('div');
            empty.style.cssText = 'color:var(--dim);font-size:11px;font-style:italic';
            empty.textContent = '(no bits exposed by server)';
            c.appendChild(empty);
            return;
        }
        bits.forEach(d => {
            const lbl = document.createElement('label');
            const cb  = document.createElement('input');
            cb.type = 'checkbox';
            cb.dataset.bit = d.name;
            cb.checked = !!(set && set.has(d.name));
            cb.onchange = syncBitMutualExclusion;
            lbl.appendChild(cb);
            lbl.appendChild(document.createTextNode(d.label || d.name));
            c.appendChild(lbl);
        });
    }

    // A bit can be in Accept OR Reject but not both — disable the twin
    // checkbox when its sibling is checked, so the user can't accidentally
    // submit a contradictory filter.  Called on every checkbox change AND
    // after the lists are built (to seed the disabled state from saved
    // filter values).
    function syncBitMutualExclusion(){
        const accept = $('cut-accept-list');
        const reject = $('cut-reject-list');
        if (!accept || !reject) return;
        const accBits = new Set(
            Array.from(accept.querySelectorAll('input[type="checkbox"]'))
                .filter(cb => cb.checked).map(cb => cb.dataset.bit));
        const rejBits = new Set(
            Array.from(reject.querySelectorAll('input[type="checkbox"]'))
                .filter(cb => cb.checked).map(cb => cb.dataset.bit));
        const setDisabled = (cb, otherSet) => {
            const taken = otherSet.has(cb.dataset.bit);
            cb.disabled = taken;
            const lbl = cb.parentElement;
            if (lbl) lbl.classList.toggle('disabled', taken);
        };
        accept.querySelectorAll('input[type="checkbox"]')
            .forEach(cb => setDisabled(cb, rejBits));
        reject.querySelectorAll('input[type="checkbox"]')
            .forEach(cb => setDisabled(cb, accBits));
    }

    // Populate every form field from a {time?, integral?, height?, quality_bits?}
    // shape.  Used by openCutDialog (current runtime filter) and the Reset
    // button (file-config defaults).
    function populateForm(filter){
        const f = filter || {};
        const fillAxis = (axis, idMin, idMax) => {
            const r = f[axis] || {};
            $(idMin).value = r.min != null ? r.min : '';
            $(idMax).value = r.max != null ? r.max : '';
        };
        fillAxis('time',     'cut-time-min',     'cut-time-max');
        fillAxis('integral', 'cut-integral-min', 'cut-integral-max');
        fillAxis('height',   'cut-height-min',   'cut-height-max');

        const qb = f.quality_bits || {};
        buildBitList('cut-accept-list', new Set(qb.accept || []));
        buildBitList('cut-reject-list', new Set(qb.reject || []));
        syncBitMutualExclusion();
    }

    function openCutDialog(){
        setDialogOpen('cut', true);
        $('cut-status-msg').textContent = '';
        populateForm(histConfig.waveform_filter);
    }

    function readAxis(idMin, idMax){
        const a = $(idMin).value, b = $(idMax).value;
        const out = {};
        if (a !== '') out.min = parseFloat(a);
        if (b !== '') out.max = parseFloat(b);
        return Object.keys(out).length ? out : null;
    }

    function readBitNames(containerId){
        return Array.from($(containerId).querySelectorAll('input[type="checkbox"]'))
            .filter(cb => cb.checked).map(cb => cb.dataset.bit);
    }

    function buildFilter(){
        const f = {};
        const t = readAxis('cut-time-min',     'cut-time-max');     if (t) f.time     = t;
        const i = readAxis('cut-integral-min', 'cut-integral-max'); if (i) f.integral = i;
        const h = readAxis('cut-height-min',   'cut-height-max');   if (h) f.height   = h;
        const acc = readBitNames('cut-accept-list');
        const rej = readBitNames('cut-reject-list');
        if (acc.length || rej.length) f.quality_bits = {accept: acc, reject: rej};
        return f;
    }

    // Local redraw — pulls overlays from histConfig.waveform_filter and
    // the cut-show state.  Server isn't touched.  Bypasses the
    // histogram refresh throttle so toggles and saves feel responsive —
    // without this, showHistograms() may early-return for ~1s and the
    // overlays don't update.
    function redrawAll(){
        lastHistModule = '';
        if (selectedModule) showWaveform(selectedModule);
        else redrawGeo();
    }

    // Server roundtrip + force redraw once histConfig is refreshed.
    // Used by both the Save button and the apply-toggle so changes
    // appear immediately even if the user hasn't toggled show.
    function refreshAfterServer(){
        fetchConfigAndApply().then(redrawAll, redrawAll);
    }

    $('cut-settings-btn').onclick = openCutDialog;
    const closeCutDialog = wireDialogClose('cut', 'cut-cancel');

    // Reset → restore the file-config filter from monitor_config.json
    // (snapshotted server-side as `waveform_filter_default`).  The user
    // still has to click Save to commit; this just repopulates the form.
    $('cut-reset').onclick = () => populateForm(histConfig.waveform_filter_default);

    $('cut-apply-btn').onclick = () => {
        const body = {
            waveform_filter:        buildFilter(),
            waveform_filter_active: $('cut-apply').checked
        };
        $('cut-status-msg').textContent = 'Saving…';
        postJson('/api/hist_config', body).then(d => {
            if (d.error) {
                $('cut-status-msg').textContent = 'Error: ' + d.error;
                return;
            }
            closeCutDialog();
            refreshAfterServer();
        }).catch(() => {
            $('cut-status-msg').textContent = 'Request failed';
        });
    };

    // "apply" toggle: immediate server POST (flips peak_filter.enable).
    $('cut-apply').onchange = function(){
        postJson('/api/hist_config', {waveform_filter_active: this.checked})
            .then(refreshAfterServer).catch(() => {});
    };

    // "show" toggle: client-side overlay only — no server roundtrip.
    $('cut-show').onchange = redrawAll;
}
