// Online mode: WebSocket + ring buffer
function updateRingSelector() {
    fetch('/api/ring').then(r => r.json()).then(data => {
        const sel = document.getElementById('ring-select');
        const prev = sel.value;
        sel.innerHTML = '';
        const ring = data.ring || [];
        for (let i = ring.length - 1; i >= 0; i--) {
            const o = document.createElement('option');
            o.value = ring[i];
            o.textContent = `Sample ${ring[i]}` + (i === ring.length - 1 ? ' (latest)' : '');
            sel.appendChild(o);
        }
        if (!autoFollow && prev && ring.includes(parseInt(prev))) sel.value = prev;
        else if (ring.length) sel.value = ring[ring.length - 1];
    });
}

function setEtStatus(connected, waiting, retries) {
    const el = document.getElementById('et-status');
    if (connected) {
        el.textContent = '● Connected';
        el.style.color = THEME.success;
    } else if (waiting) {
        el.textContent = `● Waiting for ET (${retries||'...'})`;
        el.style.color = THEME.warn;
    } else {
        el.textContent = '● Disconnected';
        el.style.color = THEME.danger;
    }
}

function updateFollowStatus() {
    const el = document.getElementById('follow-status');
    if (!el) return;
    if (autoFollow) {
        el.style.display = 'none';
    } else {
        el.textContent = `⏸ Paused at ${sampleLabel()} — click or press F to resume`;
        el.style.display = '';
    }
}

// MONITOR STATUS — poll server for DAQ livetime + beam (server shells out
// to caget for each).  Thresholds + poll interval + units come from the
// server config in applyConfig().  monitorStatusEnabled is true when at
// least one of livetime / beam-energy / beam-current has a command set.
let monitorStatusPollMs=5000;
let livetimeHealthy=90, livetimeWarning=80;
let livetimeEnabled=false, livetimeUnit='%';
let beamEnergyEnabled=false, beamEnergyUnit='MeV';
let beamCurrentEnabled=false, beamCurrentUnit='nA';
let beamCurrentTripWarn=null;     // null = no warn threshold configured
let monitorStatusEnabled=false;
let monitorStatusTimer=null;
function pollMonitorStatus(){
    fetch('/api/monitor_status').then(r=>r.json()).then(d=>{
        const ltEl=document.getElementById('livetime-display');
        if(ltEl){
            const lt=d.livetime||{};
            const ts  =(lt.ts>=0)      ? lt.ts       : null;
            const meas=(lt.measured>=0)? lt.measured : null;
            if(ts==null && meas==null){
                if(livetimeEnabled){
                    ltEl.style.display='';
                    ltEl.textContent='DAQ Livetime: N/A';
                    ltEl.style.color=THEME.textDim;
                } else {
                    ltEl.style.display='none';
                }
            } else {
                ltEl.style.display='';
                const parts=[];
                if(ts!=null)   parts.push(ts.toFixed(1)+livetimeUnit+' (TS)');
                if(meas!=null) parts.push(meas.toFixed(1)+livetimeUnit+' (DSC)');
                ltEl.textContent='DAQ Livetime: '+parts.join(' / ');
                // Color by the worse of the two so a sick channel still flags red.
                const worst=Math.min(ts??meas, meas??ts);
                ltEl.style.color=worst>=livetimeHealthy?THEME.success
                              :worst>=livetimeWarning?THEME.warn:THEME.danger;
            }
        }

        const beamEl=document.getElementById('beam-display');
        if(beamEl){
            const beam=d.beam||{};
            const eVal=(beam.energy>=0) ? beam.energy : null;
            const iVal=(beam.current>=0)? beam.current: null;
            const showE=beamEnergyEnabled, showI=beamCurrentEnabled;
            if(!showE && !showI){
                beamEl.style.display='none';
            } else {
                beamEl.style.display='';
                const ePart=showE ? (eVal!=null ? `${eVal.toFixed(1)} ${beamEnergyUnit}` : `N/A ${beamEnergyUnit}`) : '';
                const iPart=showI ? (iVal!=null ? `${iVal.toFixed(2)} ${beamCurrentUnit}` : `N/A ${beamCurrentUnit}`) : '';
                let txt='Beam: ';
                if(showE && showI)      txt += `${ePart} @ ${iPart}`;
                else if(showE)          txt += ePart;
                else                    txt += iPart;
                beamEl.textContent=txt;
                const tripped = (showI && beamCurrentTripWarn!=null
                                 && iVal!=null && iVal < beamCurrentTripWarn);
                const missing = (showE && eVal==null) || (showI && iVal==null);
                beamEl.style.color = tripped ? THEME.danger
                                   : missing ? THEME.textDim
                                             : THEME.success;
            }
        }
    }).catch(()=>{});
}
function startMonitorStatusPolling(){
    if(monitorStatusTimer) return;
    if(!monitorStatusEnabled) return;     // nothing configured server-side
    pollMonitorStatus();
    monitorStatusTimer=setInterval(pollMonitorStatus,monitorStatusPollMs);
}
function stopMonitorStatusPolling(){
    if(monitorStatusTimer){clearInterval(monitorStatusTimer);monitorStatusTimer=null;}
    for(const id of ['livetime-display','beam-display']){
        const el=document.getElementById(id);
        if(el) el.style.display='none';
    }
}

function connectWebSocket() {
    const proto = location.protocol === 'https:' ? 'wss:' : 'ws:';
    ws = new WebSocket(`${proto}//${location.host}`);

    ws.onopen = () => {
        // Advertise auto-report support; only clients that send this are
        // in the server's dispatchCapture candidate pool.
        try {
            ws.send(JSON.stringify({
                type: 'client_hello',
                capabilities: ['auto_report'],
            }));
        } catch(e) {}
    };
    ws.onclose = () => {
        setTimeout(connectWebSocket, 2000);
    };
    ws.onmessage = (evt) => {
        try {
            const msg = JSON.parse(evt.data);
            const now = Date.now();
            if (msg.type === 'new_event') {
                setEtStatus(true);  // receiving events means ET is connected
                if (autoFollow && now - lastEventFetch > refreshEventMs) {
                    lastEventFetch = now;
                    loadLatestEvent();
                }
                if (now - lastRingFetch > refreshRingMs) {
                    lastRingFetch = now;
                    updateRingSelector();
                }
                if (now - lastOccFetch > refreshHistMs) {
                    lastOccFetch = now;
                    if(histEnabled) { fetchOccupancy(); fetchClHist(); }
                    if(activeTab==='physics') fetchPhysics();
                    if(activeTab==='gem') fetchGemAccum();
                    if(activeTab==='cluster') fetchGemResiduals();
                }
                // gem_apv tab is per-event; loadEventData (called by
                // loadLatestEvent above) hooks into activeTab and refetches
                // once currentEvent has been updated.
            } else if (msg.type === 'status') {
                setEtStatus(msg.connected, msg.waiting, msg.retries);
            } else if (msg.type === 'hist_cleared') {
                occData={}; occTcutData={}; occTotal=0;
                resetClusterHists();
                lastHistModule = '';   // bypass refresh throttle
                // Use showWaveform (not showHistograms) so the waveform
                // plot's cut-range shapes also refresh when the peak
                // filter changes.  showWaveform calls showHistograms +
                // redrawGeo internally; cached samples avoid a re-fetch.
                if (selectedModule) showWaveform(selectedModule);
                else                redrawGeo();
                clearPhysicsFrontend();
                if(activeTab==='gem') fetchGemAccum();
            } else if (msg.type === 'hist_config_updated') {
                // Server's peak_filter changed (any client could have edited).
                // Pull fresh config so histConfig.waveform_filter and the
                // "apply" checkbox stay in sync, then redraw histograms and
                // the geo (color metric uses the time filter).
                fetchConfigAndApply();
                if (selectedModule) showHistograms(selectedModule);
                redrawGeo();
            } else if (msg.type === 'lms_event') {
                if (now - lastLmsFetch > refreshLmsMs) {
                    lastLmsFetch = now;
                    if(activeTab==='lms'){ fetchLmsSummary(); refreshSelectedLmsHistory(); }
                }
            } else if (msg.type === 'lms_cleared') {
                lmsSummaryData=null; resetLmsSelection();
                if(activeTab==='lms'){ geoLms(); updateLmsTable(); }
            } else if (msg.type === 'epics_event') {
                if (now - lastEpicsFetch > refreshEpicsMs) {
                    lastEpicsFetch = now;
                    if(activeTab==='epics'){
                        fetchEpicsChannels();
                        fetchEpicsLatest();
                        fetchAllEpicsSlots();
                    }
                }
            } else if (msg.type === 'epics_cleared') {
                clearEpicsFrontend();
            } else if (msg.type === 'mode_changed') {
                if (msg.mode && msg.mode !== mode) {
                    clearFrontend();
                    fetchConfigAndApply();
                }
            } else if (msg.type === 'autoclear_done') {
                // Server-side autoclear wiped histograms / lms / epics; reset
                // the whole UI so every connected tab clears together.  The
                // preceding *_cleared handlers already ran; clearFrontend is
                // idempotent on top of them.
                clearFrontend();
            } else if (msg.type === 'capture_request') {
                // Server picked this client as the on-demand reporter for
                // the run; only the selected connection receives this.
                handleCaptureRequest(msg);
            } else if (msg.type === 'auto_capture_done') {
                // Authoritative end-of-flow signal: drop the reporter badge.
                autoSetReporting(false);
                const sb = document.getElementById('status-bar');
                if (sb) sb.textContent = msg.posted
                    ? `Auto report posted (run ${msg.run||'?'})`
                    : `Auto report saved locally (run ${msg.run||'?'})`;
            } else if (msg.type === 'gem_threshold_updated') {
                // Another viewer (or this one) changed the GEM σ — sync this
                // tab's input; hits[] update with the next event.
                syncGemApvZsSigmaInput(msg.zs_sigma);
                if (gemApvCalib) gemApvCalib.zs_sigma = msg.zs_sigma;
            } else if (msg.type === 'gem_apv_full_event') {
                // Server captured a new "monitoring event" — one where the
                // firmware bypassed online ZS, so every channel of every
                // APV was read out.  Only the gem_apv tab in 'Latest full-
                // readout' source cares; gemApvOnLiveEvent gates by source
                // and pause state.
                if (activeTab === 'gem_apv') gemApvOnLiveEvent(msg.seq || 0, 'full_event');
            }
        } catch (e) {}
    };
}

