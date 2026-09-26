// viewer.js — Orchestrator: shared state, init, tab switching, event navigation

// ── State ─────────────────────────────────────────────────────────────
let modules=[], totalEvents=0, currentEvent=1;
let currentEventNumber=0, currentTriggerBits=0;  // DAQ event number + trigger from last loaded event
let currentEventKind='physics';                  // 'physics' | 'sync' | 'epics' | 'prestart' | ...
let currentRunNumber=0;                          // DAQ run number from last loaded event (0 if unknown)
let triggerBitsDef=[];  // [{bit, mask, name, label}, ...]
let triggerTypeDef=[];  // [{type, tag, name, label, primary_bit}, ...]
// per-tab trigger filter masks: { tabName: {accept: mask, reject: mask} }
const tabTrigFilter={};
function trigFilter(){ return tabTrigFilter[activeTab]||(tabTrigFilter[activeTab]={accept:0,reject:0}); }
let filterActive=false, filteredIndices=null, filteredCount=0;
let eventChannels={};
let selectedModule=null, hoveredModule=null;
const PC=['#00b4d8','#ff6b6b','#51cf66','#ffd43b','#cc5de8','#ff922b','#20c997','#f06595'];
const CRATE_NAME={0x80:'adchycal1',0x82:'adchycal2',0x84:'adchycal3',0x86:'adchycal4',0x88:'adchycal5',0x8a:'adchycal6',0x8c:'adchycal7',
    0x01:'PRadTS',0x04:'PRadROC_1',0x05:'PRadROC_2',0x06:'PRadROC_3',0x07:'PRadSRS_1',0x08:'PRadSRS_2'};
function crateName(r){return CRATE_NAME[r]||`ROC 0x${r.toString(16)}`;}
let histEnabled=false, histConfig={};
let mode='idle';    // 'idle', 'file', or 'online'
let etAvailable=false, fileAvailable=false;
let sourceCaps={type:'',has_waveforms:true,has_peaks:true,has_clusters:true,has_epics:true};
let ws=null;        // WebSocket connection (always connected)
let autoFollow=true; // auto-load latest event
let lastEventFetch=0, lastHistFetch=0, lastRingFetch=0, lastOccFetch=0, lastLmsFetch=0;
let refreshEventMs=200, refreshRingMs=500, refreshHistMs=2000, refreshLmsMs=2000;

// occupancy data (fetched once per file load when histograms enabled)
let occData={}, occTcutData={}, occTotal=0;

let activeTab='dq';  // key of TAB_PANELS

// ── Plotly shared config — PL is a GETTER so it tracks the active theme ──
Object.defineProperty(window, 'PL', { get: () => plotlyLayout() });
const PC2={responsive:true,displayModeBar:false};
const PC_EPICS={responsive:true,displayModeBar:true,
    modeBarButtonsToRemove:['sendDataToCloud','lasso2d','select2d'],
    displaylogo:false};

// POST `body` as JSON; resolves to the parsed JSON reply.
function postJson(url, body){
    return fetch(url,{method:'POST',headers:{'Content-Type':'application/json'},
        body:JSON.stringify(body)}).then(r=>r.json());
}

// ── Plotly plot registry ──────────────────────────────────────────────
// All Plotly divs register here with their tab and empty-state title.
const plotRegistry=[];  // [{id, tab, title, config}]

function registerPlot(id, tab, title, config){
    plotRegistry.push({id, tab, title, config: config||PC2});
}

// Built on every call so the empty state follows the active theme.
function emptyPlotLayout(p){
    const layout = plotlyLayout();
    if(p.title) layout.title={text:p.title, font:{size:10,color:THEME.textMuted}};
    return layout;
}

function initRegisteredPlots(){
    for(const p of plotRegistry)
        Plotly.newPlot(p.id, [], emptyPlotLayout(p), p.config);
}

// Blank a registered plot back to its empty state.
function resetPlot(id){
    const p=plotRegistry.find(q=>q.id===id);
    if(p) Plotly.react(id, [], emptyPlotLayout(p), p.config);
}

function resizePlots(...ids){
    for(const id of ids) try{Plotly.Plots.resize(id);}catch(e){}
}

function resizePlotsForTab(tab){
    resizePlots(...plotRegistry.filter(p=>p.tab===tab).map(p=>p.id));
}

function resizeAllPlots(){
    resizePlots(...plotRegistry.map(p=>p.id));
}

function redrawGeo(){
    if(activeTab==='cluster') geoCluster();
    else if(activeTab==='lms') geoLms();
    else geoDq();
}

function geoHandleClick(cx,cy){
    const m=hitTest(cx,cy);
    if(!m){
        // click on empty canvas — deselect
        if(activeTab==='cluster'){
            selectedCluster=-1;
            document.getElementById('cl-select').value='all';
            geoCluster(); updateClusterTable(); showClusterDetail();
        } else if(activeTab==='lms'){
            resetLmsSelection();
            updateLmsTable(); geoLms();
        } else {
            selectedModule=null;
            document.getElementById('detail-header').innerHTML=
                '<div class="empty-msg">Click a module to view details</div>';
            resetDqPlots('');
            geoDq();
        }
        return;
    }
    if(activeTab==='cluster'){
        selectedModule=null;
        const found=clusterOfModule(modules.indexOf(m));
        selectedCluster=(found<0||selectedCluster===found)?-1:found;
        document.getElementById('cl-select').value=selectedCluster>=0?selectedCluster:'all';
        geoCluster(); updateClusterTable(); showClusterDetail();
    } else if(activeTab==='lms'){
        const idx=modules.indexOf(m);
        lmsSelectedModule=idx;
        fetchLmsHistory(idx, m.n);
        updateLmsTable();
        geoLms();
    } else {
        showWaveform(m);
    }
}

// ── Event loading (file and online mode) ──────────────────────────────
let eventRequestId = 0;  // increments on each fetch, stale responses ignored

// Build sample label: "Sample 100 (Evt. 99)"
function sampleLabel(){
    const evn=currentEventNumber?` (Evt. ${currentEventNumber})`:'';
    return `Sample ${currentEvent}${evn}`;
}

function decodeTriggerBits(bits){
    if(!bits) return '';
    const names=[];
    for(const d of triggerBitsDef){
        if(bits & (1<<d.bit)) names.push(d.name);
    }
    let s=` trig=0x${bits.toString(16)}`;
    if(names.length) s+=` [${names.join('+')}]`;
    return s;
}

// Human-readable label for non-physics event types.  Sync/Epics/control
// samples are kept in the index so the EPICS tab + control bookkeeping see
// them, but they decode to empty FADC data — show the kind instead of "0
// channels, no trigger".
function eventKindLabel(kind){
    switch(kind){
        case 'sync':     return 'Sync event (TI sync, no FADC readout)';
        case 'epics':    return 'EPICS slow-control event';
        case 'prestart': return 'Prestart (run start marker)';
        case 'go':       return 'Go (run go marker)';
        case 'end':      return 'End (run end marker)';
        case 'control':  return 'Control event';
        case 'unknown':  return 'Unknown event type';
        default:         return '';
    }
}

// Header-centre run-number readout. Hidden until we get a non-zero
// run_number from /api/event/*; updates whenever we advance to a new
// event with a fresh run number.
function updateRunDisplay(){
    const el = document.getElementById('run-display');
    if(!el) return;
    if(currentRunNumber > 0){
        el.textContent = 'Run #'+currentRunNumber;
        el.style.display = '';
    } else {
        el.style.display = 'none';
    }
}

function updateStatusBar(){
    const modeTag = mode === 'online' ? ' [LIVE]' : '';
    const trig = decodeTriggerBits(currentTriggerBits);
    const label = sampleLabel();
    const sb = document.getElementById('status-bar');

    if(currentEventKind && currentEventKind !== 'physics'){
        // Non-physics sample — replace the per-tab readout summary with the
        // event kind so the user knows it's not a missed/empty readout.
        sb.textContent = `${label}: ${eventKindLabel(currentEventKind)}${modeTag}`;
        return;
    }

    if(activeTab==='cluster'){
        const nc=clusterData?clusterData.clusters?clusterData.clusters.length:0:0;
        const nh=clusterData?clusterData.hits?Object.keys(clusterData.hits).length:0:0;
        sb.textContent=`${label}: ${nc} clusters, ${nh} hit modules${trig}${modeTag}`;
    } else if(activeTab==='lms'){
        const lmsN=lmsSummaryData?lmsSummaryData.events||0:0;
        sb.textContent=`${label} | LMS: ${lmsN} events${modeTag}`;
    } else {
        const nch = Object.keys(eventChannels).length;
        const npk = Object.values(eventChannels).reduce((s,c) => s + (c.pk||[]).length, 0);
        sb.textContent=`${label}: ${nch} channels, ${npk} peaks${trig}${modeTag}`;
    }
}

// ── Trigger filter — accept/reject events by trigger bits ─────────────
// Each trigger bit has 3 states: unchecked (ignore), accept (green), reject (red).
// Accept: event must have at least one accepted bit set.
// Reject: event must NOT have any rejected bit set.
// If no accept bits selected, accept-all (only reject mask applies).

// Set a trigger checkbox to state st (0 ignore, 1 accept, 2 reject) and
// colour its label to match.
function setTrigCbState(cb, st){
    cb.dataset.state=String(st);
    cb.checked=st===1;
    cb.indeterminate=st===2;
    cb.parentElement.className=['','trig-accept','trig-reject'][st];
}

function buildTriggerFilterUI(){
    const bar=document.getElementById('trigger-filter-bar');
    const container=document.getElementById('trigger-filter-checks');
    if(!container) return;
    container.innerHTML='';
    if(!triggerBitsDef.length){ bar.style.display='none'; return; }
    bar.style.display='flex';

    for(const d of triggerBitsDef){
        const bit = d.bit;
        if(bit===undefined) continue;
        const mask = 1 << bit;
        const name = d.label || d.name;
        const lbl=document.createElement('label');
        lbl.title=`${name} — FP bit ${bit}`;
        const cb=document.createElement('input');
        cb.type='checkbox';
        cb.dataset.bit=bit;
        cb.dataset.mask=mask;
        cb.dataset.state='0';
        cb.indeterminate=false;
        cb.checked=false;
        cb.addEventListener('click', ()=>{
            // state is read by saveTrigFilterToTab right away; the visuals
            // wait until the browser has applied its own checkbox toggle
            const st=(parseInt(cb.dataset.state)+1)%3;
            cb.dataset.state=String(st);
            setTimeout(()=>setTrigCbState(cb,st),0);
            saveTrigFilterToTab();
        });
        lbl.appendChild(cb);
        lbl.appendChild(document.createTextNode(`${name}(${bit})`));
        container.appendChild(lbl);
    }

    document.getElementById('trig-filter-clear').onclick=()=>{
        for(const cb of container.querySelectorAll('input[type="checkbox"]'))
            setTrigCbState(cb,0);
        saveTrigFilterToTab();
    };
}

function saveTrigFilterToTab(){
    const tf=trigFilter();
    tf.accept=0; tf.reject=0;
    for(const cb of document.querySelectorAll('#trigger-filter-checks input[type="checkbox"]')){
        const m=parseInt(cb.dataset.mask);
        const st=cb.dataset.state;
        if(st==='1') tf.accept|=m;
        else if(st==='2') tf.reject|=m;
    }
    if(mode==='file' && totalEvents>0) loadEvent(currentEvent);
}

function restoreTrigFilterFromTab(){
    const tf=trigFilter();
    for(const cb of document.querySelectorAll('#trigger-filter-checks input[type="checkbox"]')){
        const m=parseInt(cb.dataset.mask);
        setTrigCbState(cb, tf.accept&m?1:tf.reject&m?2:0);
    }
}

function passesTriggerFilter(triggerBits){
    const tf=trigFilter();
    if(tf.reject && (triggerBits & tf.reject)) return false;
    if(tf.accept && !(triggerBits & tf.accept)) return false;
    return true;
}

// ── Auto Report mode ──────────────────────────────────────────────────
// Fully server-driven.  The server detects run boundaries (END /
// PRESTART control events, run-number flip on physics events as a
// fallback), picks one alive WS client, and sends a 'capture_request'.
// The chosen client takes screenshots and POSTs them to /api/elog/post;
// the server saves locally + optionally uploads to elog.  When the
// server's capture-aware autoclear fires, every client gets an
// autoclear_done broadcast that resets local UI in lockstep.
//
// This client only:
//   1. Reflects auto_post_enabled in the header status pill.
//   2. Lights the same pill green ("Auto-reporting…") while we're the
//      one running captures.

let autoPostEnabled=false;     // server-controlled, from /api/config
let autoIsReporting=false;     // true while we're handling a capture_request
// Mirrored from /api/config -> auto_report.partial_threshold_events.
// report.js reads it for the empty-state header banner.  0 disables the
// marker.  The default matches the server's so an unconfigured deployment
// still flags zeros.
let autoReportPartialThreshold=1000;

function autoStatusEl(){ return document.getElementById('auto-status'); }

function autoUpdateStatus(){
    const el=autoStatusEl(); if(!el) return;
    if(autoIsReporting){
        el.classList.add('reporting');
        el.classList.remove('on','off');
        el.textContent='Auto-reporting…';
        el.title='This browser is capturing + uploading the auto-report';
        return;
    }
    el.classList.remove('reporting');
    if(autoPostEnabled){
        el.classList.add('on');  el.classList.remove('off');
        el.textContent='Auto: ON';
        el.title='Auto-report ON — server will pick one connected client per run boundary';
    } else {
        el.classList.add('off'); el.classList.remove('on');
        el.textContent='Auto: OFF';
        el.title='Auto-report disabled in monitor_config.json';
    }
}

function autoSetReporting(on){
    autoIsReporting = !!on;
    autoUpdateStatus();
}

// Called from initReport() once /api/config has arrived.
function applyAutoReportConfig(cfg){
    if(!cfg) return;
    autoPostEnabled = !!cfg.enabled;
    autoUpdateStatus();
    if(typeof cfg.partial_threshold_events === 'number')
        autoReportPartialThreshold = cfg.partial_threshold_events;
}

// Manual Clear All — operator-driven, immediate.  Run-boundary
// (PRESTART / END / run-change) wipes are handled entirely server-side
// via scheduleAutoClear; this path is intentionally not gated on
// pending_capture_ so a button press always takes effect right away.
function doClearAll(){
    return Promise.all([
        fetch('/api/hist/clear').then(r=>r.json()),
        fetch('/api/lms/clear').then(r=>r.json()),
        fetch('/api/epics/clear').then(r=>r.json()),
    ]).then(clearFrontend).catch(()=>{
        document.getElementById('status-bar').textContent='Error clearing data';
    });
}

let navDirection=1;  // +1=forward, -1=backward (for trigger filter auto-skip)

// `manual` distinguishes user-initiated event swaps (file-mode arrow keys /
// prev-next / jump-to, online ring-buffer nav) from the auto-refresh path
// (online WS new_event → loadLatestEvent).  Currently only the gem_apv tab
// cares — its pause gate bypasses on manual=true.
function loadEventData(reqId, data, manual) {
    if (reqId !== eventRequestId) return;  // stale response, discard
    if (data.error) {
        document.getElementById('status-bar').textContent = data.error;
        if (mode === 'online') updateRingSelector();
        return;
    }

    // trigger filter: if event doesn't pass, auto-skip in nav direction
    const tf=trigFilter();
    if (tf.accept || tf.reject) {
        const tb = data.trigger_bits || 0;
        if (!passesTriggerFilter(tb)) {
            if (mode==='file' && !stepEvent(navDirection))
                document.getElementById('status-bar').textContent =
                    `No matching event (trigger filter active)`;
            // in online mode, just discard this event silently
            return;
        }
    }

    currentEvent = data.event;
    currentEventNumber = data.event_number || 0;
    currentTriggerBits = data.trigger_bits || 0;
    currentEventKind = data.event_kind || 'physics';
    if (data.run_number) {
        currentRunNumber = data.run_number;
        updateRunDisplay();
    }
    eventChannels = data.channels || {};
    if(mode==='online') sampleCount++;
    updateStatusBar();
    updateHeaderStats();
    if(activeTab==='cluster'){
        clusterEvent=-1; // invalidate cache
        loadClusterData(currentEvent);
    } else if(activeTab==='lms'){
        // LMS geo doesn't change per event — no redraw needed
    } else if(activeTab==='gem'){
        // Right panel is per-event-independent (efficiency cards + last-good
        // snapshot), refreshed on the histogram cadence by fetchGemAccum.
    } else if(activeTab==='gem_apv'){
        // Pause + source gate; user navigation (manual) bypasses pause.
        gemApvOnLiveEvent(currentEvent, 'event', !!manual);
    } else {
        geoDq();
    }
    if (activeTab==='dq' && selectedModule) showWaveform(selectedModule);
    updateGeoTooltip();

    // refresh histograms from server (on-demand accumulation in non-preprocessed mode)
    if(mode==='file' && !histEnabled){
        fetchOccupancy(); fetchClHist(); fetchGemResiduals();
    }
}

function loadEvent(evnum) {
    currentEvent = evnum;
    const reqId = ++eventRequestId;
    if (mode === 'file') {
        if(filteredIndices){
            const pos=filteredIndices.indexOf(evnum);
            document.getElementById('ev-input').value = pos>=0 ? pos+1 : evnum;
        } else {
            document.getElementById('ev-input').value = evnum;
        }
    }
    document.getElementById('status-bar').textContent = `Loading sample ${evnum}...`;
    fetch(`/api/event/${evnum}`).then(r => r.json()).then(d => loadEventData(reqId, d, true))
        .catch(err => { document.getElementById('status-bar').textContent = `Error: ${err}`; });
}

// Load the previous (dir=-1) or next (dir=+1) event, staying within the
// event-filter result when one is applied.  Returns false at either end.
function stepEvent(dir){
    let t;
    if(filteredIndices){
        const pos=filteredIndices.indexOf(currentEvent);
        if(pos<0) return false;
        t=filteredIndices[pos+dir];
    } else t=currentEvent+dir;
    if(t===undefined||t<1||t>totalEvents) return false;
    navDirection=dir;
    loadEvent(t);
    return true;
}

function loadLatestEvent() {
    const reqId = ++eventRequestId;
    fetch('/api/event/latest').then(r => r.json()).then(d => loadEventData(reqId, d, false))
        .catch(err => { document.getElementById('status-bar').textContent = `Error: ${err}`; });
}

// ── Draggable dividers ────────────────────────────────────────────────
function setupDivider(divId, axis, getTarget, getContainer, getOffset, minA, minB, onResize){
    const div=document.getElementById(divId);
    let active=false;
    div.addEventListener('mousedown',e=>{
        active=true; div.classList.add('active');
        document.body.style.cursor=axis==='x'?'col-resize':'row-resize';
        document.body.style.userSelect='none'; e.preventDefault();
    });
    document.addEventListener('mousemove',e=>{
        if(!active)return;
        const container=getContainer(), rect=container.getBoundingClientRect();
        const pos=axis==='x'?e.clientX-rect.left-getOffset():e.clientY-rect.top-getOffset();
        const max_=(axis==='x'?rect.width:rect.height)-getOffset()-minB;
        const val=Math.max(minA,Math.min(max_,pos));
        const target=getTarget();
        target.style.flex='none';  // override flex:1 so explicit size takes effect
        if(axis==='x') target.style.width=val+'px'; else target.style.height=val+'px';
        onResize();
    });
    document.addEventListener('mouseup',()=>{
        if(!active)return; active=false; div.classList.remove('active');
        document.body.style.cursor=''; document.body.style.userSelect='';
    });
}

// ── Tab switching ─────────────────────────────────────────────────────
// Panels shown for each tab (keys match the data-tab values); every other
// id in ALL_TAB_PANEL_IDS is hidden.  report.js also uses the table to
// scope its screenshot capture and to hide the other panels in the clone.
const TAB_PANELS={
    dq:      ['geo-panel','div-main','geo-toolbar-dq','detail-panel'],
    lms:     ['geo-panel','div-main','geo-toolbar-lms','lms-panel'],
    cluster: ['geo-panel','div-main','geo-toolbar-cl','cluster-panel'],
    gem:     ['gem-outer'],
    gem_apv: ['gem-apv-outer'],
    epics:   ['epics-outer'],
    physics: ['physics-outer'],
};
const ALL_TAB_PANEL_IDS=[...new Set(Object.values(TAB_PANELS).flat())];

function switchTab(tab, opts){
    if(tab===activeTab) return;
    activeTab=tab;
    selectedModule=null;
    restoreTrigFilterFromTab();
    // clear notification dot on the tab being opened
    if(tab==='lms') document.getElementById('lms-dot').className='tab-dot';
    if(tab==='epics') document.getElementById('epics-dot').className='tab-dot';
    document.querySelectorAll('.tab').forEach(t=>{
        t.classList.toggle('active', t.dataset.tab===tab);
    });
    const shown=TAB_PANELS[tab]||[];
    for(const id of ALL_TAB_PANEL_IDS)
        document.getElementById(id).style.display=
            !shown.includes(id)?'none':(id==='geo-panel'||id==='div-main')?'':'flex';

    // --- per-tab actions: fetch data + resize after layout settles ---
    const tabActions = {
        dq:      { hasGeo: true },
        cluster: { hasGeo: true,
                   fetch(){ loadClusterData(currentEvent); },
                   after(){ plotClHist(); plotClStatHists(); } },
        lms:     { hasGeo: true,
                   fetch(){ fetchLmsSummary(); } },
        epics:   { fetch(){ fetchEpicsChannels(); fetchEpicsLatest(); fetchAllEpicsSlots(); } },
        physics: { fetch(){ fetchPhysics(); } },
        gem:     { fetch(){ fetchGemAccum(); },
                   after(){ resizePlots(...GEM_OCC_IDS, 'gem-eff-zhist'); } },
        gem_apv: { fetch(){
                       // Honour pause across tab switches — if the operator
                       // froze the panel on event N, switching away and back
                       // shouldn't clobber that frame.  Otherwise pull from
                       // the active source (current event or latest full).
                       if (gemApvPaused && gemApvData) { renderGemApvPanels(); return; }
                       refreshGemApv(currentEvent);
                   },
                   after(){ resizeGemApv(); } },
    };
    const action = tabActions[tab] || tabActions.dq;

    if (action.fetch && !(opts && opts.skipFetch)) action.fetch();

    // after layout settles: resize geo (for geo tabs) + registered Plotly plots + custom after()
    setTimeout(()=>{
        if (action.hasGeo) resizeGeo();
        resizePlotsForTab(tab);
        if (action.after) action.after();
    }, 50);
    updateStatusBar();
}

// ── Modal dialogs: #<name>-dialog + #<name>-backdrop, shown with class 'open'
function setDialogOpen(name, on){
    for(const s of ['-backdrop','-dialog'])
        document.getElementById(name+s).classList.toggle('open', on);
}

// Close the dialog from its backdrop, its header × (#<name>-dialog-close)
// and any extra buttons; returns the close function.
function wireDialogClose(name, ...extraIds){
    const close=()=>setDialogOpen(name, false);
    for(const id of [name+'-backdrop', name+'-dialog-close', ...extraIds])
        document.getElementById(id).onclick=close;
    return close;
}

// ── Init ──────────────────────────────────────────────────────────────
function init(){
    drawColorBar(); initGeo();
    for(const id of COLORBAR_IDS) document.getElementById(id).onclick=cyclePalette;
    registerPlot('waveform-div', 'dq', null);
    registerPlot('heighthist-div','dq', 'Height Histogram');
    registerPlot('inthist-div',  'dq', 'Integral Histogram');
    registerPlot('poshist-div',  'dq', 'Position Histogram');

    // --- copy data buttons ---
    function setupCopyBtn(btnId, getData) {
        document.getElementById(btnId).onclick=()=>{
            const d=getData();
            if(!d) return;
            const text=`x: [${d.x.join(', ')}]\ny: [${d.y.join(', ')}]`;
            navigator.clipboard.writeText(text).then(()=>{
                const btn=document.getElementById(btnId);
                btn.textContent='✓'; setTimeout(()=>{btn.textContent='copy';},1000);
            });
        };
    }
    setupCopyBtn('btn-copy-wf', ()=>currentWaveform);
    setupCopyBtn('btn-copy-heighthist', ()=>currentHist['heighthist-div']);
    setupCopyBtn('btn-copy-inthist', ()=>currentHist['inthist-div']);
    setupCopyBtn('btn-copy-poshist', ()=>currentHist['poshist-div']);

    // --- dividers ---
    // 1. main vertical: geo ↔ detail
    setupDivider('div-main','x',
        ()=>document.getElementById('geo-panel'),
        ()=>document.querySelector('.main'),
        ()=>0, 300, 350, ()=>{resizeGeo();resizeAllPlots();});
    // 2-4. horizontal dividers between the 4 DQ cells
    // Each divider resizes the cell above it; offset = top of that cell relative to container
    function dqDivider(divId, cellId, prevCellIds){
        setupDivider(divId, 'y',
            ()=>document.getElementById(cellId),
            ()=>document.getElementById('detail-panel'),
            ()=>{
                // offset = header + all cells/dividers above this cell
                let off=document.getElementById('detail-header').offsetHeight;
                for(const id of prevCellIds){
                    const el=document.getElementById(id);
                    if(el) off+=el.offsetHeight;
                }
                return off;
            },
            60, 60, resizeAllPlots);
    }
    dqDivider('div-dq-1','dq-cell-heighthist',[]);
    dqDivider('div-dq-2','dq-cell-histrow',['dq-cell-heighthist','div-dq-1']);
    dqDivider('div-dq-3','dq-cell-waveform',['dq-cell-heighthist','div-dq-1','dq-cell-histrow','div-dq-2']);

    // --- tab switching ---
    document.querySelectorAll('.tab').forEach(t=>{
        t.onclick=()=>switchTab(t.dataset.tab);
    });

    // --- cluster controls ---
    document.getElementById('cl-select').onchange=e=>{
        selectedCluster=e.target.value==='all'?-1:parseInt(e.target.value);
        geoCluster(); updateClusterTable(); showClusterDetail();
    };
    document.getElementById('cl-log-scale').onchange=()=>{ if(activeTab==='cluster') geoCluster(); };

    registerPlot('cl-energy-hist',  'cluster', 'Cluster Energy');
    registerPlot('cl-rawe-hist',   'cluster', 'Raw Energy Sum');
    registerPlot('cl-nclust-hist', 'cluster', 'Clusters per Event');
    registerPlot('cl-nblocks-hist','cluster', 'Blocks per Cluster');
    for (let d = 0; d < 4; d++)
        registerPlot('gem-resid-' + d, 'cluster', null);
    for (let d = 0; d < 4; d++)
        registerPlot('gem-eff-grid-' + d, 'gem', null);
    registerPlot('gem-eff-zy', 'gem', null);
    setupCopyBtn('btn-copy-cl-hist', ()=>currentClHist);
    setupCopyBtn('btn-copy-cl-rawe', ()=>currentRawEnergyHist);
    setupCopyBtn('btn-copy-nclust', ()=>currentNclustHist);
    setupCopyBtn('btn-copy-nblocks', ()=>currentNblocksHist);

    // cluster energy / raw energy column divider (sibling of cl-stat-row).
    setupDivider('div-cl-eh','x',
        ()=>document.querySelector('.cl-hist-cell'),
        ()=>document.querySelector('.cl-hist-row'),
        ()=>0, 80, 80, ()=>resizePlots('cl-energy-hist','cl-rawe-hist'));

    // cluster stat row column divider
    setupDivider('div-cl-stat','x',
        ()=>document.querySelector('.cl-stat-cell'),
        ()=>document.querySelector('.cl-stat-row'),
        ()=>0, 80, 80, ()=>resizePlots('cl-nclust-hist','cl-nblocks-hist'));

    // waveform Stack / DAQ toggles (DAQ = firmware Mode 1/2/3 emulation: TET,
    // NSB/NSA, Vp, T markers per the FADC250 manual) — reset to match JS
    // state (browser may restore old form state)
    setWfMode(false,false);
    // Stack and DAQ are mutually exclusive — DAQ annotations don't compose
    // with overlaid events.
    document.getElementById('wf-stack').onchange=e=>{
        setWfMode(e.target.checked,false);
        if(selectedModule) showWaveform(selectedModule);
    };
    document.getElementById('btn-wf-stack-reset').onclick=()=>{
        wfStackTraces=[]; wfStackModKey='';
        if(selectedModule) showWaveform(selectedModule);
    };
    document.getElementById('wf-daq').onchange=e=>{
        setWfMode(false,e.target.checked);
        if(selectedModule) showWaveform(selectedModule);
    };

    // histogram log-scale toggles
    document.getElementById('heighthist-logy').onchange=()=>{ if(selectedModule) showHistograms(selectedModule); };
    document.getElementById('inthist-logy').onchange=()=>{ if(selectedModule) showHistograms(selectedModule); };
    document.getElementById('clhist-logy').onchange=plotClHist;
    document.getElementById('clrawe-logy').onchange=plotRawEnergyHist;
    setupCopyBtn('btn-copy-lms', ()=>currentLmsData);

    // cluster panel divider: histogram ↔ table
    setupDivider('div-cl-ht','y',
        ()=>document.getElementById('cl-hist-panel'),
        ()=>document.getElementById('cluster-panel'),
        ()=>0,
        80, 80, ()=>resizePlots('cl-energy-hist','cl-rawe-hist'));

    registerPlot('lms-plot', 'lms', 'LMS History');
    registerPlot('physics-plot',       'physics', null, PC_EPICS);
    registerPlot('moller-xy-plot',     'physics', null, PC_EPICS);
    registerPlot('hycal-xy-plot',      'physics', null, PC_EPICS);
    for(let i=0;i<EPICS_NUM_SLOTS;i++)
        registerPlot('epics-plot-'+i, 'epics', null, PC_EPICS);

    initRegisteredPlots();

    // GEM APV tab toolbar (Process / Signal Only / Shared Y / sample mask).
    setupGemApvControls();

    setupDivider('div-lms-ht','y',
        ()=>document.getElementById('lms-plot-panel'),
        ()=>document.getElementById('lms-panel'),
        ()=>0,
        80, 80, ()=>resizePlots('lms-plot'));
    document.getElementById('lms-color-metric').onchange=geoLms;
    document.getElementById('lms-log-scale').onchange=geoLms;
    document.getElementById('lms-ref-select').onchange=e=>{
        g_lmsRefIndex=parseInt(e.target.value);
        fetchLmsSummary();
        refreshSelectedLmsHistory();
    };

    // --- file mode nav ---
    document.getElementById('btn-prev').onclick=()=>stepEvent(-1);
    document.getElementById('btn-next').onclick=()=>stepEvent(1);
    document.getElementById('ev-input').onchange=e=>{
        const v=parseInt(e.target.value);
        if(filteredIndices){
            if(v>=1&&v<=filteredIndices.length) loadEvent(filteredIndices[v-1]);
        } else {
            if(v>=1&&v<=totalEvents) loadEvent(v);
        }
    };
    document.getElementById('color-metric').onchange=()=>{syncDqRange();geoDq();};
    document.getElementById('log-scale').onchange=geoDq;
    // Note: cut-show / cut-apply / cut-settings-btn handlers live in cut_dialog.js

    // --- file browser ---
    document.getElementById('btn-open').onclick = openFileDialog;
    wireDialogClose('file');
    document.getElementById('file-filter').oninput = e => filterFileList(e.target.value);
    document.addEventListener('keydown', e => {
        if (e.key === 'Escape')
            document.querySelectorAll('.file-dialog.open').forEach(d => setDialogOpen(d.id.replace(/-dialog$/, ''), false));
    });

    // --- range editing ---
    function setupRangeEdit(btnId, editId, showId, getVal, setVal, onApply) {
        const btn=document.getElementById(btnId);
        const edit=document.getElementById(editId);
        const show=document.getElementById(showId);
        let editing=false;

        function startEdit() {
            editing=true; btn.classList.add('editing'); btn.textContent='✓';
            edit.classList.add('active'); show.style.display='none';
            edit.value=getVal()||'';
            edit.focus(); edit.select();
        }
        function applyEdit() {
            if(!editing) return;
            editing=false; btn.classList.remove('editing'); btn.textContent='✎';
            edit.classList.remove('active'); show.style.display='';
            const v=parseFloat(edit.value);
            setVal(isNaN(v)?null:v);
            onApply();
        }

        btn.addEventListener('mousedown', e => e.preventDefault());
        btn.onclick=()=>{ if(editing) applyEdit(); else startEdit(); };
        edit.addEventListener('keydown',e=>{
            if(e.key==='Enter') applyEdit();
            if(e.key==='Escape'){editing=false;btn.classList.remove('editing');btn.textContent='✎';edit.classList.remove('active');show.style.display='';}
        });
        edit.addEventListener('blur',()=>{ applyEdit(); });
    }
    // Min/max editors (ids <prefix>range-{min,max}-{btn,edit,show}) of the
    // colour range stored per (tab, metric) in geoRangeOverrides.
    function bindRangeEditors(prefix, tab, metricFn, onApply){
        ['min','max'].forEach((w,i)=>setupRangeEdit(
            `${prefix}range-${w}-btn`, `${prefix}range-${w}-edit`, `${prefix}range-${w}-show`,
            ()=>getGeoRange(tab,metricFn())[i],
            v=>{ const r=getGeoRange(tab,metricFn()).slice(); r[i]=v; setGeoRange(tab,metricFn(),r[0],r[1]); },
            onApply));
    }
    bindRangeEditors('','dq',()=>document.getElementById('color-metric').value,()=>{syncDqRange();geoDq();});
    bindRangeEditors('cl-','cluster',()=>'energy',geoCluster);
    bindRangeEditors('lms-','lms',()=>document.getElementById('lms-color-metric').value,geoLms);

    // --- online mode nav ---
    document.getElementById('ring-select').onchange=e=>{
        autoFollow=false; updateFollowStatus();
        loadEvent(parseInt(e.target.value));
    };
    document.getElementById('ring-select').onfocus=()=>{ updateRingSelector(); };
    document.getElementById('follow-status').onclick=()=>{ autoFollow=true; updateFollowStatus(); loadLatestEvent(); };

    // Theme toggle — cycles dark → light → classic → dark
    const themeBtn = document.getElementById('btn-theme');
    const THEME_GLYPH = { dark: '☀', light: '☾', classic: '◉' };
    function updateThemeBtn(){
        const t = currentTheme();
        themeBtn.textContent = THEME_GLYPH[t] || '◐';
        themeBtn.title = `Theme: ${t} — click to cycle`;
    }
    updateThemeBtn();
    themeBtn.onclick = () => toggleTheme();
    onThemeChange(() => {
        updateThemeBtn();
        // Per-tab listeners (cluster/lms/gem/gem_apv/physics/epics) do
        // full Plotly.react replays so trace and title colors track the
        // active palette.  This relayout is the safety net for plots
        // whose listener early-returned (e.g. lms-plot before the user
        // clicks a module) — chrome at minimum stays in sync.
        for(const p of plotRegistry) plotlyRelayout(p.id);
        // Canvases read THEME at draw time, so force a full repaint.
        // drawColorBar covers every per-tab colorbar canvas; redrawGeo
        // hits whichever shared-geo tab is active.
        drawColorBar();
        redrawGeo();
        // DQ waveform + histograms — Plotly shapes (cut-range overlays,
        // ref lines) bake THEME.cutShade at draw time, so a relayout
        // chrome patch isn't enough.  Re-run showWaveform to regenerate
        // the layout (calls wfLayout + showHistograms internally).
        if (selectedModule) { lastHistModule = ''; showWaveform(selectedModule); }
    });

    // Clear All — resets all tabs' data for new run
    document.getElementById('btn-clear-all').onclick=()=>doClearAll();

    // Auto-report status pill — server picks the chosen client per END /
    // run-change; pill reflects auto_report.enabled and lights green
    // while we're the chosen reporter.
    autoUpdateStatus();

    // mode toggle button — opens ET dialog when going online
    document.getElementById('btn-mode-toggle').onclick=()=>{
        if(mode==='online'){
            fetch('/api/mode/file',{method:'POST'}).then(()=>fetchConfigAndApply());
        } else {
            openEtDialog();
        }
    };

    initFilterDialog();
    initCutDialog();


    // ET connect dialog
    const closeEtDialog=wireDialogClose('et','et-cancel');
    document.getElementById('et-connect').onclick=()=>{
        const cfg={
            host:    document.getElementById('et-input-host').value,
            port:    parseInt(document.getElementById('et-input-port').value)||11111,
            et_file: document.getElementById('et-input-file').value,
            station: document.getElementById('et-input-station').value,
        };
        document.getElementById('et-status-msg').textContent='Connecting...';
        postJson('/api/mode/online', cfg).then(d=>{
            if(d.error){
                document.getElementById('et-status-msg').textContent='Error: '+d.error;
            } else {
                closeEtDialog();
                fetchConfigAndApply();
            }
        }).catch(()=>{
            document.getElementById('et-status-msg').textContent='Connection failed';
        });
    };

    // geo mouse
    const tip=document.getElementById('geo-tooltip');
    function tooltipText(m){
        let t=`${m.n}  (${m.t==='G'?'PbGlass':'PbWO₄'})\n${crateName(m.roc)}  slot ${m.sl}  ch ${m.ch}`;
        if(activeTab==='lms' && lmsSummaryData){
            const idx=modules.indexOf(m);
            const md=lmsSummaryData.modules?lmsSummaryData.modules[String(idx)]:null;
            if(md){
                const rmsPct=md.mean>0?(md.rms/md.mean*100).toFixed(1):'--';
                t+=`\nLMS Mean: ${md.mean.toFixed(1)}  RMS: ${md.rms.toFixed(2)}  (${rmsPct}%)`;
                t+=`\n${md.count} pts  ${lmsState(md)!=='ok'?'⚠ WARNING':'OK'}`;
            } else { t+='\nNo LMS data'; }
        } else if(activeTab==='cluster' && clusterData){
            const idx=modules.indexOf(m);
            const energy=clusterData.hits?clusterData.hits[String(idx)]:0;
            if(energy) t+=`\nEnergy: ${energy.toFixed(1)} MeV`;
            const ci=clusterOfModule(idx);
            if(ci>=0){
                const cl=clusterData.clusters[ci];
                t+=`\nCluster #${ci} (${cl.center}, ${cl.energy.toFixed(0)} MeV)`;
            }
        } else {
            const key=`${m.roc}_${m.sl}_${m.ch}`;
            const d=eventChannels[key];
            if(d&&d.pk&&d.pk.length){
                const pks=peaksInCut(d.pk);
                const bp=tallest(pks);
                const tc=isTimeCut();
                if(bp) t+=`\nPed ${d.pm.toFixed(1)}  H ${bp.h.toFixed(0)}  Int ${bp.i.toFixed(0)}  T ${bp.t.toFixed(0)}ns  Pk ${pks.length}${tc?' (tcut)':''}`;
                else t+=`\nPed ${d.pm.toFixed(1)}  (no peaks${tc?' in time cut':''})`;
            }
            else if(d)t+=`\nPed ${d.pm.toFixed(1)}  (no peaks)`;
            const pct=occPct(key);
            if(pct!==null) t+=`\nOcc ${pct.toFixed(1)}%  (${occTotal} evts${isTimeCut()?' tcut':''})`;
            else if(histEnabled===false){
                t+=`\nOcc: not computed (enable histograms)`;
            }
        }
        return t;
    }
    // update tooltip content for currently hovered module (called on data refresh)
    updateGeoTooltip=()=>{
        if(hoveredModule) tip.textContent=tooltipText(hoveredModule);
    };
    geoCanvas.addEventListener('mousemove',e=>{
        const r=geoCanvas.getBoundingClientRect(),m=hitTest(e.clientX-r.left,e.clientY-r.top);
        if(m!==hoveredModule){
            hoveredModule=m;
            renderGeoOutlines(_geoOutlineFn, _geoDecorateFn);  // outlines only — fills unchanged
        }
        if(m){
            tip.textContent=tooltipText(m);tip.style.display='block';
            tip.style.left=(e.clientX-r.left+14)+'px';tip.style.top=(e.clientY-r.top-8)+'px';
        }else tip.style.display='none';
    });
    // clicks go through geoHandleClick (geo.js mouseup below the drag threshold)
    geoCanvas.addEventListener('mouseleave',()=>{
        hoveredModule=null;tip.style.display='none';
        renderGeoOutlines(_geoOutlineFn, _geoDecorateFn);
    });
    document.addEventListener('keydown',e=>{
        if(e.target.tagName==='INPUT'||e.target.tagName==='SELECT')return;
        if(mode==='file'){
            if(e.key==='ArrowLeft') stepEvent(-1);
            if(e.key==='ArrowRight') stepEvent(1);
        } else {
            if(e.key==='ArrowLeft'||e.key==='ArrowRight'){
                // navigate ring buffer
                const sel=document.getElementById('ring-select');
                const opts=[...sel.options].map(o=>parseInt(o.value));
                const cur=parseInt(sel.value);
                const idx=opts.indexOf(cur);
                if(e.key==='ArrowRight'&&idx>0){sel.value=opts[idx-1];autoFollow=false;updateFollowStatus();loadEvent(opts[idx-1]);}
                if(e.key==='ArrowLeft'&&idx<opts.length-1){sel.value=opts[idx+1];autoFollow=false;updateFollowStatus();loadEvent(opts[idx+1]);}
            }
            if(e.key==='f'||e.key==='F'){autoFollow=true;updateFollowStatus();loadLatestEvent();}
        }
    });

    // always connect WebSocket (for mode_changed and clear notifications)
    connectWebSocket();

    // load config and init mode
    fetchConfigAndApply();
}
window.addEventListener('DOMContentLoaded',init);

