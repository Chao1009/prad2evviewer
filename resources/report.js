// report.js — Report generation: markdown body + per-tab PNG screenshots.
//
// One screenshot per tab (Waveform / LMS / Clustering / GEM / APV / EPICS / Physics)
// instead of capturing individual plots — captureTabScreenshot() switches to
// the tab, swaps every visible <canvas> and Plotly plot with a freshly
// rendered <img>, then composites the cloned DOM through an SVG
// <foreignObject> wrapper so toolbars, tables, and labels round-trip with
// their CSS intact.
//
// Driven by the server's capture_request message (handleCaptureRequest), which
// generates the report on the chosen client and posts it to the elog.
//
// Depends on globals from viewer.js (accessed at runtime, not load time).

// ── Registry ──────────────────────────────────────────────────────────
const reportRegistry=[];
let elogConfig={url:'',logbook:'',author:'',tags:[]};
let reportAttachments=[];  // [{data (base64), filename, caption, type}]

function registerReportSection(section){
    reportRegistry.push(section);
    reportRegistry.sort((a,b)=>a.order-b.order);
}

function addAttachment(dataUrl,filename,caption){
    if(!dataUrl) return;
    const b64=dataUrl.split(',')[1];
    if(b64) reportAttachments.push({data:b64,filename,caption,type:'image/png'});
}

// ── Whole-tab screenshot ──────────────────────────────────────────────
// Strategy: switch to the target tab, pre-render every <canvas> + Plotly
// plot to a data URL on the live DOM (cloneNode loses canvas pixels and
// Plotly's internal SVG state), clone the panel, replace those nodes with
// <img> elements, and ship the result through an SVG <foreignObject>.  All
// same-origin stylesheets are inlined in a <style> tag so layout + theme
// survive the trip.
const TAB_SETTLE_MS = 1000;   // post-switch settling for layout + data fetch
const THEME_SETTLE_MS = 250;  // settling for async Plotly re-renders

// Tabs included in the report, in order. Keys must match data-tab values.
const REPORT_TABS = [
    {tab:'dq',      title:'Waveform Data',     filename:'tab_dq.png',      prep:_prepDqTab},
    {tab:'lms',     title:'Gain Monitoring',   filename:'tab_lms.png',     prep:_prepLmsTab},
    {tab:'cluster', title:'Clustering',        filename:'tab_cluster.png', prep:_prepClusterTab},
    {tab:'gem',     title:'GEM Detectors',     filename:'tab_gem.png',     prep:null},
    {tab:'gem_apv', title:'GEM APV Waveforms', filename:'tab_gem_apv.png', prep:null},
    {tab:'epics',   title:'EPICS Slow Control',filename:'tab_epics.png',   prep:_prepEpicsTab},
    {tab:'physics', title:'Physics',           filename:'tab_physics.png', prep:_prepPhysicsTab},
];

function _wait(ms){ return new Promise(r=>setTimeout(r,ms)); }

// Force the DQ tab to a deterministic capture state and return a restore
// function.  Five things move:
//   1) color-metric  -> 'occupancy'   (so the geo legend matches the title)
//   2) focus-pbwo4   -> unchecked     (deterministic Lead-Glass dim state)
//   3) selectedModule -> null         (no per-channel detail panel leak)
//   4) detail-header / four detail Plotly divs / peaks tables -> placeholder
//   5) repaint via geoDq() so the swap lands before the screenshot
// clearFrontend() deliberately keeps selectedModule across server-side
// autoclear so the live UI stays responsive — that's why we have to
// wipe-and-restore here per-capture instead of relying on the server's
// hist_cleared broadcast.
function _prepDqTab(){
    const sel           = document.getElementById('color-metric');
    const focus         = document.getElementById('focus-pbwo4');
    const detailHeader  = document.getElementById('detail-header');
    const peaksTbody    = document.getElementById('peaks-tbody');
    const peaksTbodyDaq = document.getElementById('peaks-tbody-daq');

    const prevMetric    = sel ? sel.value : '';
    const prevFocus     = focus ? focus.checked : false;
    const prevFocusVar  = geoFocusPbWO4;
    const prevSelected  = selectedModule;
    const prevWf        = currentWaveform;
    const prevHist      = currentHist;
    const prevHeader    = detailHeader ? detailHeader.innerHTML : '';
    const prevPeaks     = peaksTbody ? peaksTbody.innerHTML : '';
    const prevPeaksDaq  = peaksTbodyDaq ? peaksTbodyDaq.innerHTML : '';

    if(sel && prevMetric!=='occupancy'){
        sel.value='occupancy';
        syncDqRange();
    }
    if(focus && prevFocus){
        focus.checked=false;
    }
    geoFocusPbWO4=false; selectedModule=null;
    if(detailHeader)
        detailHeader.innerHTML =
            '<div class="empty-msg">No module selected (auto-report capture)</div>';
    try{ resetDqPlots(''); }catch(e){ /* missing plotly divs on a stripped page — non-fatal */ }
    geoDq();

    return ()=>{
        if(sel && sel.value!==prevMetric){
            sel.value=prevMetric;
            syncDqRange();
        }
        if(focus && focus.checked!==prevFocus){
            focus.checked=prevFocus;
        }
        geoFocusPbWO4=prevFocusVar; selectedModule=prevSelected;
        currentWaveform=prevWf; currentHist=prevHist;
        if(detailHeader)  detailHeader.innerHTML=prevHeader;
        if(peaksTbody)    peaksTbody.innerHTML=prevPeaks;
        if(peaksTbodyDaq) peaksTbodyDaq.innerHTML=prevPeaksDaq;
        // If the user had a focused channel, re-render its waveform +
        // histograms + peaks table from the live event stream so the live
        // UI returns to the prior state once the capture finishes.
        if(prevSelected){
            try{ showWaveform(prevSelected); }catch(e){}
        }
        geoDq();
    };
}

// Force the LMS tab to RMS/Mean (warning view) for the screenshot. Reference
// is left to whatever lmsSummaryData was loaded with — _summaryLms() runs
// before captureTabScreenshot and pulls that data with LMS3 as the
// reference.  We only flip the dropdown so the toolbar in the screenshot
// reads "LMS3" alongside the rms_frac coloring.  Restore on return.
function _prepLmsTab(){
    const metric=document.getElementById('lms-color-metric');
    const refSel=document.getElementById('lms-ref-select');
    const prevMetric=metric.value;
    const prevRef=refSel.value;
    if(prevMetric!=='rms_frac'){
        metric.value='rms_frac';
        geoLms();
    }
    const lms3=_lmsRefOptionValue();
    if(lms3 && refSel.value!==lms3){
        refSel.value=lms3;
        geoLms();
    }
    return ()=>{
        if(metric.value!==prevMetric){
            metric.value=prevMetric;
            geoLms();
        }
        if(refSel.value!==prevRef){
            refSel.value=prevRef;
            geoLms();
        }
    };
}

// Value of the 'LMS3' option in the LMS reference dropdown, or '' if absent.
function _lmsRefOptionValue(){
    for(const o of document.getElementById('lms-ref-select').options)
        if(o.textContent==='LMS3') return o.value;
    return '';
}

// Re-run the tab's plot* calls once it is visible.  refreshDataForReport ran
// them while the panel was display:none, so Plotly laid them out against a
// zero-size container (autorange / ticks / legend / heatmap aspect differ
// from the live monitor), and Plotly.Plots.resize only updates dimensions
// without re-running the trace pipeline.  Nothing to restore, but the no-op
// restore makes captureTabScreenshot wait THEME_SETTLE_MS for the re-render.
// Shared by the Cluster, EPICS and Physics preps.
function _prepClusterTab(){
    // Includes the shared geo (HyCal occupancy on the left) so it
    // matches the right-half plots' freshly re-react'd traces.
    try{
        plotClHist();
        plotRawEnergyHist();
        plotClStatHists();
        plotGemResiduals();
        geoCluster();
    }catch(e){ /* fall through — captureTabScreenshot still snapshots what it has */ }
    return ()=>{};
}

function _prepEpicsTab(){
    for(let i=0;i<EPICS_NUM_SLOTS;i++) plotEpicsSlot(i);
    return ()=>{};
}

function _prepPhysicsTab(){
    plotEnergyAngle(); plotMollerXY(); plotHycalXY();
    return ()=>{};
}

function _gatherCss(){
    const out=[];
    for(const sh of document.styleSheets){
        try{ for(const r of sh.cssRules) out.push(r.cssText); }
        catch(e){ /* CORS-protected sheet — skip */ }
    }
    return out.join('\n');
}

// Mirror live form-control state into attributes so XMLSerializer captures
// it (cloneNode copies attributes, not properties).
function _freezeFormState(root){
    for(const inp of root.querySelectorAll('input')){
        if(inp.type==='checkbox'||inp.type==='radio'){
            if(inp.checked) inp.setAttribute('checked','');
            else inp.removeAttribute('checked');
        }else{
            inp.setAttribute('value', inp.value);
        }
    }
    for(const sel of root.querySelectorAll('select')){
        for(const opt of sel.options){
            if(opt.selected) opt.setAttribute('selected','');
            else opt.removeAttribute('selected');
        }
    }
}

function _pathFromRoot(el, root){
    const path=[];
    let n=el;
    while(n && n!==root){
        const p=n.parentNode;
        if(!p) return null;
        path.unshift([...p.children].indexOf(n));
        n=p;
    }
    return n===root ? path : null;
}

function _resolveByPath(root, path){
    let n=root;
    for(const i of path){
        if(!n||!n.children||i<0||i>=n.children.length) return null;
        n=n.children[i];
    }
    return n;
}

// Capture an SVG-foreignObject snapshot of the .main panel for the given tab.
// Returns a PNG data URL, or null on failure.  Uses the active theme — no
// theme flipping — so the screenshot matches what the user sees on screen.
async function captureTabScreenshot(tab){
    const prevTab = activeTab;
    if(tab !== activeTab){
        // skipFetch: refreshDataForReport already loaded every tab's data and
        // ran the corresponding plot* call. Letting switchTab re-fire its
        // action.fetch races the screenshot — the second response can land
        // mid-render, leaving the captured PNG with empty plots.
        switchTab(tab, {skipFetch:true});
        await _wait(TAB_SETTLE_MS);
    } else {
        // Same-tab path: the prior section's prep (e.g. _prepLmsTab calling
        // geoLms) may have repainted the SHARED geo canvas with the wrong
        // tab's coloring.  switchTab's setTimeout-based redraw doesn't run
        // here because tab===activeTab short-circuits, so we'd snapshot the
        // stale pixels.  Force a redraw matching the current tab; canvas
        // 2D draws are synchronous so toDataURL would already pick up the
        // new bitmap, but we yield once anyway in case redrawGeo schedules
        // any deferred work (e.g. a Plotly resize on a co-mounted plot).
        redrawGeo();
        await _wait(0);
    }

    const tabSpec = REPORT_TABS.find(t=>t.tab===tab);
    const restore = tabSpec && tabSpec.prep ? tabSpec.prep() : null;
    if(restore) await _wait(THEME_SETTLE_MS);

    let dataUrl=null;
    let svgUrl=null;
    try{
        const root = document.querySelector('.main');
        const rect = root.getBoundingClientRect();
        const W = Math.max(800, Math.ceil(rect.width));
        const H = Math.max(600, Math.ceil(rect.height));

        // 1) Snapshot every visible <canvas> + Plotly plot to data URLs.
        // Restrict the search to the active tab's panel(s) so a stale
        // canvas inside a hidden panel can't sneak in if its bounding
        // rect happens to be non-zero (Plotly + display:none can leave
        // last-rendered dimensions on the .js-plotly-plot wrapper).
        const wantedPanels = TAB_PANELS[tab] || [];
        const scopes = wantedPanels.map(id=>document.getElementById(id))
                                   .filter(el=>el);
        const scopeRoots = scopes.length ? scopes : [root];
        const replacements=[];  // [{path, src, w, h}]
        const seen=new Set();
        for(const scope of scopeRoots){
            for(const c of scope.querySelectorAll('canvas')){
                if(seen.has(c)) continue; seen.add(c);
                const r=c.getBoundingClientRect();
                if(r.width===0||r.height===0) continue;
                const path=_pathFromRoot(c, root);
                if(!path) continue;
                try{ replacements.push({path, src:c.toDataURL('image/png'), w:r.width, h:r.height}); }
                catch(e){ /* tainted canvas */ }
            }
            for(const p of scope.querySelectorAll('.js-plotly-plot')){
                if(seen.has(p)) continue; seen.add(p);
                const r=p.getBoundingClientRect();
                if(r.width===0||r.height===0) continue;
                const path=_pathFromRoot(p, root);
                if(!path) continue;
                try{
                    const url=await Plotly.toImage(p,{format:'png',width:r.width,height:r.height});
                    replacements.push({path, src:url, w:r.width, h:r.height});
                }catch(e){ /* uninitialised plot */ }
            }
        }

        // 2) Clone DOM, freeze form state, swap snapshotted nodes for <img>.
        const clone = root.cloneNode(true);
        _freezeFormState(clone);
        // Belt-and-suspenders: force-hide every non-active tab panel in
        // the clone. switchTab already sets the live DOM correctly, but
        // a stray re-show (theme replay, Plotly straggler, etc.) would
        // otherwise leak the wrong outer into the cloned screenshot.
        for(const panelId of ALL_TAB_PANEL_IDS){
            if(wantedPanels.includes(panelId)) continue;
            const el=clone.querySelector('#'+panelId);
            if(el) el.style.display='none';
        }
        clone.style.width=W+'px';
        clone.style.height=H+'px';
        clone.style.background=THEME.bg;
        for(const rep of replacements){
            const twin=_resolveByPath(clone, rep.path);
            if(!twin || !twin.parentNode) continue;
            const img=document.createElement('img');
            img.setAttribute('src', rep.src);
            const w=Math.round(rep.w), h=Math.round(rep.h);
            img.setAttribute('width', String(w));
            img.setAttribute('height', String(h));
            img.setAttribute('style',
                `display:block;width:${w}px;height:${h}px;border:none;margin:0;padding:0`);
            twin.parentNode.replaceChild(img, twin);
        }

        // 3) Inline stylesheets + override transient overlays.
        const cssText=_gatherCss();
        const overrides=
            `body{background:${THEME.bg}!important;}`+
            `.main{background:${THEME.bg}!important;overflow:hidden;}`+
            `.geo-tooltip,.backdrop,.file-dialog,.progress-overlay{display:none!important;}`;

        const xhtml=new XMLSerializer().serializeToString(clone);
        const svg=
            `<svg xmlns="http://www.w3.org/2000/svg" width="${W}" height="${H}">`+
              `<foreignObject x="0" y="0" width="${W}" height="${H}">`+
                `<div xmlns="http://www.w3.org/1999/xhtml" `+
                     `style="width:${W}px;height:${H}px;background:${THEME.bg};color:${THEME.text};`+
                     `font-family:'Consolas','SF Mono',monospace;">`+
                  `<style>${cssText}\n${overrides}</style>`+
                  xhtml+
                `</div>`+
              `</foreignObject>`+
            `</svg>`;

        // 4) Render SVG → PNG via Blob URL.
        const blob=new Blob([svg], {type:'image/svg+xml'});
        svgUrl=URL.createObjectURL(blob);
        await new Promise((resolve, reject)=>{
            const img=new Image();
            img.onload=()=>{
                const out=document.createElement('canvas');
                out.width=W; out.height=H;
                const ctx=out.getContext('2d');
                ctx.fillStyle=THEME.bg;
                ctx.fillRect(0,0,W,H);
                ctx.drawImage(img, 0, 0, W, H);
                try{ dataUrl=out.toDataURL('image/png'); resolve(); }
                catch(e){ reject(e); }
            };
            img.onerror=()=>reject(new Error('SVG image load failed'));
            img.src=svgUrl;
        });
    }catch(e){
        console.error('captureTabScreenshot failed for', tab, e);
    }finally{
        if(svgUrl) URL.revokeObjectURL(svgUrl);
        if(restore) restore();
        if(prevTab !== activeTab) switchTab(prevTab, {skipFetch:true});
    }
    return dataUrl;
}

// ── Markdown table helper ─────────────────────────────────────────────
function mdTable(headers,rows,alignments){
    const aligns=alignments||headers.map(()=>'l');
    const sepMap={l:':---',r:'---:',c:':---:'};
    let md='| '+headers.join(' | ')+' |\n';
    md+='| '+aligns.map(a=>sepMap[a]||':---').join(' | ')+' |\n';
    for(const row of rows)
        md+='| '+row.join(' | ')+' |\n';
    return md+'\n';
}

// ── Per-tab text summaries ────────────────────────────────────────────
// Each section grabs whatever metadata the live globals expose and formats
// a short markdown summary; the visual is provided separately by
// captureTabScreenshot(tab).
async function _summaryDq(){
    if(occTotal>0) return `Total events: ${occTotal}\n\n`;
    // Placeholder ensures the section under the screenshot is never bare,
    // so a reviewer can't mistake "missing data" for "section intentionally
    // empty".  Same pattern repeats in the other _summary* fns below.
    return `_No waveform events accumulated (occTotal=0)._\n\n`;
}

async function _summaryLms(){
    // refresh with LMS3 ref so the warning summary is meaningful
    const lms3=parseInt(_lmsRefOptionValue());
    const refQ=lms3>=0?`?ref=${lms3}`:'';
    let d=null;
    try{ d=await fetch(`/api/lms/summary${refQ}`).then(r=>r.json()); }catch(e){}
    if(d) {
        lmsSummaryData=d;
        // refreshDataForReport's fetchLmsSummary populated the table
        // with the user's default ref; _summaryLms refetches with LMS3
        // ref so the warn table in the markdown body matches the
        // screenshot.  Re-render the DOM table too — otherwise the LMS
        // tab screenshot shows the stale "default-ref" rows (or "No
        // LMS data" if the autoclear had wiped lmsSummaryData and the
        // first fetch returned no modules yet).
        updateLmsTable();
    }
    if(!d || !d.modules || !Object.keys(d.modules).length){
        const tf=d&&d.trigger||{};
        const trigMask=`accept=0x${(tf.trigger_accept||0).toString(16)} reject=0x${(tf.trigger_reject||0).toString(16)}`;
        return `LMS events received: ${d?d.events||0:0} (trigger mask = ${trigMask})\n\n`;
    }
    const allEntries=lmsSortedEntries(d.modules);
    const driftEntries=allEntries.filter(e=>e._st==='drift');
    const warnEntries =allEntries.filter(e=>e._st==='warn');
    const okEntries   =allEntries.filter(e=>e._st==='ok');
    // Always include every drift+warn row, plus top 5 OK rows for context.
    const tableEntries=[...driftEntries, ...warnEntries, ...okEntries.slice(0,5)];

    let s=`LMS events: ${d.events||0} | Modules: ${allEntries.length} (ref: LMS3)\n\n`;
    // Per-type ranges + suppress list — show in the header so the reviewer
    // knows what counted as a flag.
    const fmtBand=(lo,hi)=> (lo!=null&&hi!=null)?`[${lo.toFixed(2)}, ${hi.toFixed(2)}]`:'?';
    const wBand=fmtBand(d.drift_low_w, d.drift_high_w);
    const gBand=fmtBand(d.drift_low_g, d.drift_high_g);
    const suppressed=(d.drift_suppress_types||[]);
    const suppressTxt=suppressed.length?` | suppressed: ${suppressed.join(', ')}`:'';
    if(driftEntries.length){
        s+=`**GAIN DRIFT**: ${driftEntries.length} module(s) outside `+
           `W ${wBand}, G ${gBand}${suppressTxt}`;
        if(d.drift_lamp_scale!=null) s+=` (lamp scale via ${d.drift_lamp_ref||'?'} = ${d.drift_lamp_scale.toFixed(3)})`;
        s+=`\n\n`;
    } else if(d.drift_enabled){
        s+=`Gain drift: 0 modules out of band — W ${wBand}, G ${gBand}${suppressTxt}`;
        if(d.drift_lamp_ref) s+=`, lamp ref ${d.drift_lamp_ref}`;
        s+=`\n\n`;
    }
    s+=`Warnings (stability): **${warnEntries.length}** / ${allEntries.length} modules\n\n`;
    if(tableEntries.length){
        const fmtState=st=>st==='drift'?'**DRIFT**':st==='warn'?'**WARN**':'OK';
        s+=mdTable(
            ['Module','Mean','RMS','RMS/Mean %','Drift','Count','Status'],
            tableEntries.map(e=>[
                e.name, e.mean.toFixed(1), e.rms.toFixed(2),
                (e.mean>0?(e.rms/e.mean*100).toFixed(1):'--')+'%',
                (e.drift!=null?e.drift.toFixed(2):'--'),
                e.count, fmtState(e._st)
            ]),['l','r','r','r','r','r','l']
        );
        if(okEntries.length>5)
            s+=`*Showing ${driftEntries.length} drift + ${warnEntries.length} warn + top 5 of ${okEntries.length} OK modules*\n\n`;
    }
    return s;
}

async function _summaryCluster(){
    if(clHistBins && clHistBins.some(b=>b>0)) {
        const sum=clHistBins.reduce((a,b)=>a+b,0);
        return `Cluster events: ${clHistEvents||0} | Histogram entries: ${sum}\n\n`;
    }
    // raw_energy may have entries even when cluster_energy doesn't (events
    // passed cluster_trigger but produced zero clusters); mention it so the
    // 'Raw Energy Sum' bin in the screenshot isn't mysterious.
    const raw=rawEnergyBins ? rawEnergyBins.reduce((a,b)=>a+b,0) : 0;
    return raw>0
        ? `_No reconstructed clusters yet (raw-energy fills = ${raw})._\n\n`
        : `_No clustering events accumulated._\n\n`;
}

async function _summaryGemApv(){
    if(!gemApvData || !gemApvData.enabled)
        return `_No GEM APV snapshot available (full-readout event not seen since last clear)._\n\n`;
    const ndets=(gemApvData.detectors||[]).length;
    const napvs=(gemApvData.apvs||[]).length;
    let s=`Event: ${gemApvCurrentEvent} | Detectors: ${ndets} | APVs: ${napvs}`;
    if(gemApvData.zs_sigma) s+=` | ZS σ: ${gemApvData.zs_sigma}`;
    return s+'\n\n';
}

async function _summaryPhysics(){
    const data=physicsData, ml=mollerData, hxy=hycalXyData;
    if((!data||!data.events) && (!ml||!ml.total_events) && (!hxy||!hxy.total_events))
        return `_No physics events accumulated (energy/angle, Møller, HyCal-XY all empty)._\n\n`;
    const evts=data?.events || ml?.total_events || hxy?.total_events || 0;
    let s=`Events: ${evts}`;
    if(data?.beam_energy) s+=` | Beam: ${data.beam_energy.toFixed(2)} MeV`;
    if(data?.hycal_z) s+=` | HyCal z: ${(data.hycal_z/1000).toFixed(2)}m`;
    if(ml) s+=` | Møller: ${ml.moller_events}`;
    if(hxy) s+=` | HyCalXY: ${hxy.events}`;
    return s+'\n\n';
}

const _SECTION_SUMMARIES = {
    dq:      _summaryDq,
    lms:     _summaryLms,
    cluster: _summaryCluster,
    gem_apv: _summaryGemApv,
    // GEM, EPICS: screenshot is enough; the EPICS channel table duplicated
    // it in text and dominated the post body, so no text summary.
    physics: _summaryPhysics,
};

// Register the per-tab sections in display order.
REPORT_TABS.forEach((t, i)=>{
    registerReportSection({
        id:t.tab, title:t.title, order:i+1,
        generate: async ()=>{
            const summary = _SECTION_SUMMARIES[t.tab];
            const summaryMd = summary ? await summary().catch(()=>'') : '';
            const img = await captureTabScreenshot(t.tab);
            if(img) addAttachment(img, t.filename, t.title);
            let md=`## ${t.title}\n\n`;
            if(summaryMd) md+=summaryMd;
            if(img) md+=`![${t.title}](${t.filename})\n\n`;
            return md;
        }
    });
});

// Pre-fetch live data for every tab so captureTabScreenshot doesn't have
// to wait for in-flight async fetches inside switchTab. Each fetch* call
// returns a promise that resolves after its plot* call has run, so by the
// time Promise.all finishes the page state matches what each tab would
// show if the user navigated to it.
async function refreshDataForReport(){
    const fetches=[fetchOccupancy(), fetchClHist(), fetchGemResiduals(), fetchGemAccum(),
        fetchLmsSummary(), fetchEpicsChannels(), fetchEpicsLatest(), fetchAllEpicsSlots(),
        fetchPhysics()];
    if(currentEvent>0) fetches.push(fetchGemApvData(currentEvent));
    await Promise.all(fetches);
    // Plotly.react schedules an async DOM update; give the browser a frame
    // to commit the latest plots before we start switching tabs.
    await _wait(THEME_SETTLE_MS);
}

// ── Report generation core ────────────────────────────────────────────
// Generate the report. Returns {md, attachments, low_data, partial} or null.
async function generateReport(reportBy,runNumber){
    if(!modules.length){
        alert('No data loaded. Please load data before generating a report.');
        return null;
    }
    const statusBar=document.getElementById('status-bar');
    const prevStatus=statusBar.textContent;
    statusBar.textContent='Generating report...';
    try{
        await refreshDataForReport();
        reportAttachments=[];
        const ts=new Date().toLocaleString();
        // online: use the server's events_processed (total field on
        // /api/occupancy → occTotal).  sampleCount is the chosen
        // client's UI auto-follow counter; on a fresh / backgrounded /
        // recently-cleared tab it can read 1-2 even when the server
        // has 10k+ events for the run, which is what the report is
        // actually summarising.
        const samples=mode==='online'?occTotal:totalEvents;
        // Empty-state detection.  refreshDataForReport has just landed
        // the three counters as globals, so the read is synchronous.
        // _suspicious -> every accumulator that has a per-tab summary
        //                line was empty at capture.
        // _partial    -> BOTH occupancy/recon samples AND cluster
        //                events are below partialThr.  Deliberately
        //                an AND, not an OR: a cluster-rich run with
        //                undersampled occupancy (or vice-versa) is
        //                still informative enough not to flag.  Set
        //                partial_threshold_events = 0 to disable.
        const clusterEvts = clHistEvents;
        const lmsEvts     = (lmsSummaryData && lmsSummaryData.events) || 0;
        const _suspicious = (samples===0 && clusterEvts===0 && lmsEvts===0);
        const partialThr  = autoReportPartialThreshold;
        const _partial   = !_suspicious && partialThr>0 &&
                            samples<partialThr && clusterEvts<partialThr;
        const runStr=runNumber?String(runNumber).padStart(6,'0'):'';
        const titleRun=runStr?`Run ${runStr}: `:'';
        let header=`# ${titleRun}PRad-II HyCal Monitor Report\n\n`;
        header+=`- **Generated:** ${ts}\n`;
        header+=`- **Samples:** ${samples}\n`;
        if(runNumber) header+=`- **DAQ Run:** ${runNumber}\n`;
        if(reportBy) header+=`- **Report by:** ${reportBy}\n`;
        if(_suspicious){
            header+=`\n> **WARNING: No data captured for this report.**\n`+
                    `> events_processed = ${samples}, cluster_events = ${clusterEvts}, lms_events = ${lmsEvts}.\n`+
                    `> Likely causes: server restarted mid-run, schedule trigger fired before the run accumulated data, or accumulators were manually cleared just before the capture.\n`+
                    `> Inspect the screenshots manually before treating this as a run state summary.\n\n`;
        } else if(_partial){
            header+=`\n> **Note:** Partial-run report (samples=${samples}, cluster events=${clusterEvts}).\n`+
                    `> The schedule trigger fired before the run had accumulated typical statistics; subsequent reports for this run will supersede this one.\n\n`;
        }
        let sectionsMd='';
        for(const entry of reportRegistry){
            try{
                const section=await entry.generate();
                if(section) sectionsMd+=section;
            }catch(err){
                sectionsMd+=`## ${entry.title}\n\n*Error: ${err.message}*\n\n`;
            }
        }
        // append LMS warn summary to header (data available after LMS section runs)
        if(lmsSummaryData&&lmsSummaryData.modules){
            const warns=Object.values(lmsSummaryData.modules)
                .filter(m=>lmsState(m)!=='ok').map(m=>m.name)
                .sort((a,b)=>{
                    const ta=a.startsWith('W')?0:1, tb=b.startsWith('W')?0:1;
                    if(ta!==tb) return ta-tb;
                    return a.localeCompare(b,undefined,{numeric:true});
                });
            if(warns.length)
                header+=`- **Gain Monitoring Warnings (${warns.length}):** ${warns.join(', ')}\n`;
            else
                header+=`- **Gain Monitoring Warnings:** None\n`;
        }
        let md=header+`\n---\n\n`+sectionsMd;
        md+=`---\n*PRad-II HyCal Online Monitor — Report generated ${ts}*\n`;
        statusBar.textContent=prevStatus;
        return {md, attachments:reportAttachments,
                low_data:_suspicious, partial:_partial};
    }catch(err){
        statusBar.textContent=`Report error: ${err.message}`;
        return null;
    }
}

// ── Elog XML helpers (used by handleCaptureRequest) ────────────────────
function escXml(s){
    return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;')
        .replace(/>/g,'&gt;').replace(/"/g,'&quot;');
}

function buildElogXml(title,logbook,author,tags,body,attachments){
    const parts=['<?xml version="1.0" encoding="UTF-8"?>','<Logentry>',
        `  <created>${new Date().toISOString()}</created>`,
        `  <Author><username>${escXml(author)}</username></Author>`,
        `  <title>${escXml(title)}</title>`,
        `  <body type="text"><![CDATA[${body}]]></body>`,
        '  <Logbooks>'];
    for(const lb of logbook.split(','))
        parts.push(`    <logbook>${escXml(lb.trim())}</logbook>`);
    parts.push('  </Logbooks>');
    if(tags&&tags.length){
        parts.push('  <Tags>');
        for(const t of tags) parts.push(`    <tag>${escXml(t.trim())}</tag>`);
        parts.push('  </Tags>');
    }
    if(attachments&&attachments.length){
        parts.push('  <Attachments>');
        for(const a of attachments)
            parts.push('    <Attachment>',
                `      <caption>${escXml(a.caption)}</caption>`,
                `      <filename>${escXml(a.filename)}</filename>`,
                `      <type>${escXml(a.type)}</type>`,
                `      <data encoding="base64">${a.data}</data>`,
                '    </Attachment>');
        parts.push('  </Attachments>');
    }
    parts.push('</Logentry>');
    return parts.join('\n');
}

// ── On-demand capture handler ─────────────────────────────────────────
// Called by online.js's WS handler when the server sends
// {type:"capture_request", request_id, run, reason}.  Only the chosen
// client receives this message — every other tab ignores it.  No client-
// side dedup or rate-limit logic; the server is the single gatekeeper.
function elogAvailable(){ return !!(elogConfig && elogConfig.url); }

function autoReportTitle(runNumber, lowData){
    const baseTitle='PRad2 Event Monitor Auto Report';
    // [NO-DATA] prefix makes a bad report unmistakable in the elog
    // titles list before anyone opens the entry.  See the _suspicious
    // flag in generateReport for the trigger condition.
    const prefix = lowData ? '[NO-DATA] ' : '';
    return runNumber ? `${prefix}Run #${runNumber}: ${baseTitle}`
                     : `${prefix}${baseTitle}`;
}

async function handleCaptureRequest(msg){
    const runNumber = msg && msg.run     || 0;
    const reason    = msg && msg.reason  || 'auto';
    const requestId = msg && msg.request_id || '';
    const sb = document.getElementById('status-bar');

    if(!elogAvailable() || !modules.length){
        if(sb) sb.textContent = `Auto-report skipped: ${
            !elogAvailable() ? 'elog not configured' : 'no modules loaded'}`;
        return;
    }

    autoSetReporting(true);
    if(sb) sb.textContent = `Auto-report (${reason}, run ${runNumber}): capturing…`;

    try{
        const reportBy = elogConfig.author || 'auto';
        const report   = await generateReport(reportBy, runNumber||'');
        if(!report) {
            if(sb) sb.textContent = 'Auto-report failed: report generation';
            return;
        }
        // generateReport sets low_data when every per-tab accumulator was
        // empty at capture time — used here for the title prefix and
        // threaded through the POST so the server filenames + summary
        // record agree with the body banner.
        const lowData   = !!report.low_data;
        const fullTitle = autoReportTitle(runNumber, lowData);
        // Tags must come from the logbook's enumerated set; the reason
        // ("end" / "run-change") is informational only, so we surface
        // it in the body rather than as a tag (it would fail PRADLOG's
        // tag-enum schema validation).
        const tags = (elogConfig.tags||[]).slice();
        const body = `*Auto-posted: ${reason}*\n\n`+
                     report.md.replace(/!\[[^\]]*\]\([^)]+\)\n*/g,'');
        const xml  = buildElogXml(fullTitle, elogConfig.logbook||'',
                                  reportBy, tags, body, report.attachments);

        if(sb) sb.textContent = `Auto-report (run ${runNumber}): uploading…`;
        const r = await postJson('/api/elog/post', {xml, auto:true, run_number:runNumber,
                                                    request_id:requestId, low_data:lowData});
        if(r.ok && r.dry_run){
            if(sb) sb.textContent =
                `Auto-report saved (dry-run): ${r.saved_xml||r.saved_dir||fullTitle}`;
        } else if(r.ok && r.skipped){
            if(sb) sb.textContent = `Auto-report skipped: ${r.detail||fullTitle}`;
        } else if(r.ok){
            if(sb) sb.textContent = `Auto-report posted: ${fullTitle}`;
        } else {
            if(sb) sb.textContent = `Auto-report failed: ${
                r.error || ('HTTP '+(r.status||'?'))}`;
        }
        // The server's auto_capture_done WS broadcast will also clear the
        // reporter badge on every client (including this one). Clear it
        // here too to keep the local UI snappy.
    } catch(e){
        if(sb) sb.textContent = `Auto-report error: ${e.message}`;
    } finally {
        autoSetReporting(false);
    }
}

// Called from config.js with the /api/config data.  Stashes the elog
// config (nested under auto_report) so handleCaptureRequest has
// author/logbook/tags to draw from, and forwards the auto_report config
// block on to viewer.js's applyAutoReportConfig.
function initReport(data){
    const ar = data && data.auto_report;
    if(ar && ar.elog && ar.elog.url) elogConfig = ar.elog;
    applyAutoReportConfig(ar);
}
