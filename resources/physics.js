// physics.js — Physics tab: HyCal cluster XY (left) + Møller XY (top right) + energy vs angle (bottom right)
//
// Depends on globals from viewer.js: PL, PC_EPICS, triggerBitsDef

let physicsData=null, mollerData=null, hycalXyData=null;

function fetchEnergyAngle(){
    return fetch('/api/physics/energy_angle').then(r=>r.json()).then(data=>{
        physicsData=data;
        plotEnergyAngle();
    }).catch(()=>{});
}

function fetchMoller(){
    return fetch('/api/physics/moller').then(r=>r.json()).then(data=>{
        mollerData=data;
        plotMollerXY();
    }).catch(()=>{});
}

function fetchHycalXY(){
    return fetch('/api/physics/hycal_xy').then(r=>r.json()).then(data=>{
        hycalXyData=data;
        plotHycalXY();
    }).catch(()=>{});
}

function fetchPhysics(){
    return Promise.all([fetchEnergyAngle(), fetchMoller(), fetchHycalXY()]);
}

// Heatmap trace of a row-major nx*ny histogram with bin centres x0+(i+0.5)*dx,
// y0+(j+0.5)*dy.  logZ plots log10 of the counts (empty bins blank); the hover
// text always shows the raw counts.
function physHeatmap(bins, nx, ny, x0, dx, y0, dy, logZ, hover, logTitle){
    const z=[], text=[];
    for(let iy=0;iy<ny;iy++){
        const row=bins.slice(iy*nx,(iy+1)*nx);
        z.push(logZ?row.map(v=>v>0?Math.log10(v):null):row);
        text.push(row.map(v=>String(v)));
    }
    return {
        z,
        x:Array.from({length:nx},(_,i)=>x0+(i+0.5)*dx),
        y:Array.from({length:ny},(_,i)=>y0+(i+0.5)*dy),
        type:'heatmap', colorscale:'Hot', reversescale:false,
        hovertemplate:hover, text,
        colorbar:{title:logZ?logTitle:'counts',titleside:'right',
            titlefont:{size:10,color:THEME.textDim},tickfont:{size:9,color:THEME.textDim}},
    };
}

// Layout of the equal-scale X/Y (mm) hit maps.
function physXYLayout(title, shapes){
    return {...PL,
        title:{text:title,font:{size:11,color:THEME.text}},
        xaxis:{...PL.xaxis,title:'X (mm)',scaleanchor:'y',scaleratio:1},
        yaxis:{...PL.yaxis,title:'Y (mm)'},
        margin:{l:50,r:70,t:30,b:35},
        shapes,
    };
}

// ep elastic scattering: E' = E / (1 + (E/Mp)*(1 - cos(theta)))
function elasticEp(beamE, thetaDeg){
    const Mp=938.272;
    const th=thetaDeg*Math.PI/180;
    return beamE/(1+(beamE/Mp)*(1-Math.cos(th)));
}

function plotEnergyAngle(){
    const div='physics-plot';
    if(!physicsData||!physicsData.bins||!physicsData.bins.length||!physicsData.nx){
        Plotly.react(div,[],{...PL,title:{text:'Energy vs Angle — No data',font:{size:12,color:THEME.textDim}}},PC_EPICS);
        document.getElementById('physics-stats').textContent='';
        return;
    }
    const d=physicsData;
    const logZ=document.getElementById('physics-logz').checked;
    const showElastic=document.getElementById('physics-elastic').checked;

    const traces=[physHeatmap(d.bins,d.nx,d.ny,d.angle_min,d.angle_step,d.energy_min,d.energy_step,logZ,'θ=%{x:.2f}° E=%{y:.0f} MeV: %{text}<extra></extra>','log₁₀(counts)')];

    if(showElastic && d.beam_energy>0){
        const ex=[],ey=[];
        for(let th=d.angle_min+0.1;th<=d.angle_max;th+=0.05){
            const e=elasticEp(d.beam_energy,th);
            if(e>=d.energy_min&&e<=d.energy_max){ex.push(th);ey.push(e);}
        }
        traces.push({x:ex,y:ey,mode:'lines',
            line:{color:THEME.success,width:2,dash:'dot'},
            name:`ep elastic (${d.beam_energy.toFixed(2)} MeV)`,
            hovertemplate:'θ=%{x:.2f}° E=%{y:.0f} MeV<extra>ep elastic</extra>'});
    }

    Plotly.react(div,traces,{...PL,
        title:{text:`Energy vs Angle (${d.events} evts)`,font:{size:12,color:THEME.text}},
        xaxis:{...PL.xaxis,title:'Scattering Angle (deg)'},
        yaxis:{...PL.yaxis,title:'Energy (MeV)'},
        margin:{l:55,r:80,t:30,b:40},
        showlegend:showElastic,
        legend:{x:0.7,y:0.95,font:{size:10,color:THEME.textDim},bgcolor:'rgba(0,0,0,0)'},
        shapes:refShapes('energy_angle'),
    },PC_EPICS);

    const ml=mollerData;
    let stats=`${d.events} evts | beam: ${d.beam_energy>0?d.beam_energy.toFixed(2):'?'} MeV`;
    if(ml) stats+=` | Møller: ${ml.moller_events}`;
    document.getElementById('physics-stats').textContent=stats;
}

function plotMollerXY(){
    const div='moller-xy-plot';
    const d=mollerData;
    if(!d||!d.xy_bins||!d.xy_bins.length||!d.xy_nx){
        Plotly.react(div,[],{...PL,title:{text:'Møller XY — No data',font:{size:12,color:THEME.textDim}}},PC_EPICS);
        return;
    }
    const logZ=document.getElementById('physics-logz').checked;
    const cuts=d.cuts||{};
    const fmtA=v=>v!=null?v.toFixed(2):'?';
    // Trigger tag: which trigger stream feeds this monitor (X17 runs take
    // Møllers from the 2-cluster trigger).  Resolve the accept mask against
    // triggerBitsDef for the display labels; fall back to server-sent names.
    // accept==0 (accept-all) shows no tag.
    const acc=(d.trigger&&d.trigger.trigger_accept)||0;
    let trigNames=[];
    if(acc && triggerBitsDef.length)
        trigNames=triggerBitsDef.filter(t=>acc&(1<<t.bit)).map(t=>t.label||t.name);
    else if(d.trigger_accept_names) trigNames=d.trigger_accept_names;
    const trigTxt=trigNames.length?`[${trigNames.join('+')}] `:'';
    const cutTxt=`${trigTxt}θ∈[${fmtA(cuts.angle_min)},${fmtA(cuts.angle_max)}]° Esum±${((cuts.energy_tolerance||0.1)*100).toFixed(0)}%`;

    // θ ring overlay: convert the moller angle window into HyCal-plane radii
    // via r = dz · tan(θ), centered at (target_x, target_y).  dz is the
    // target→HyCal lever arm — same geometry the server uses to compute the
    // per-cluster theta in app_state.cpp.
    const shapes=[...(refShapes('moller_xy')||[])];
    const target=d.target||[0,0,0];
    const dz=((d.hycal_z!=null?d.hycal_z:0)-target[2]);
    if(dz>0 && cuts.angle_min!=null && cuts.angle_max!=null){
        const cx=target[0], cy=target[1];
        const ringColor=THEME.accent;
        [cuts.angle_min,cuts.angle_max].forEach(thDeg=>{
            const r=dz*Math.tan(thDeg*Math.PI/180);
            shapes.push({
                type:'circle', xref:'x', yref:'y',
                x0:cx-r, x1:cx+r, y0:cy-r, y1:cy+r,
                line:{color:ringColor,width:1.2,dash:'dash'},
                fillcolor:'rgba(0,0,0,0)',
            });
        });
    }

    Plotly.react(div,[physHeatmap(d.xy_bins,d.xy_nx,d.xy_ny,d.xy_x_min,d.xy_x_step,d.xy_y_min,d.xy_y_step,logZ,'x=%{x:.1f} y=%{y:.1f} mm: %{text}<extra></extra>','log₁₀')],physXYLayout(`Møller XY (${d.moller_events} evts) ${cutTxt}`,shapes),PC_EPICS);
}

function plotHycalXY(){
    const div='hycal-xy-plot';
    const d=hycalXyData;
    if(!d||!d.xy_bins||!d.xy_bins.length||!d.xy_nx){
        Plotly.react(div,[],{...PL,title:{text:'HyCal Cluster Hits — No data',font:{size:12,color:THEME.textDim}}},PC_EPICS);
        return;
    }
    const logZ=document.getElementById('physics-logz').checked;
    const c=d.cuts||{};
    const fracPct=((c.energy_frac_min||0.9)*100).toFixed(0);
    const cutTxt=`Ncl=${c.n_clusters||1}, E≥${fracPct}% Eb, blocks∈[${c.nblocks_min||0},${c.nblocks_max||0}]`;

    Plotly.react(div,[physHeatmap(d.xy_bins,d.xy_nx,d.xy_ny,d.xy_x_min,d.xy_x_step,d.xy_y_min,d.xy_y_step,logZ,'x=%{x:.1f} y=%{y:.1f} mm: %{text}<extra></extra>','log₁₀')],physXYLayout(`HyCal Cluster Hits (${d.events} evts) ${cutTxt}`,refShapes('hycal_xy')),PC_EPICS);
}

function clearPhysicsFrontend(){
    physicsData=null; mollerData=null; hycalXyData=null;
    Plotly.react('physics-plot',[],{...PL},PC_EPICS);
    Plotly.react('moller-xy-plot',[],{...PL},PC_EPICS);
    Plotly.react('hycal-xy-plot',[],{...PL},PC_EPICS);
    document.getElementById('physics-stats').textContent='';
}

function initPhysics(data){
    document.getElementById('physics-logz').onchange=()=>{plotEnergyAngle();plotMollerXY();plotHycalXY();};
    document.getElementById('physics-elastic').onchange=plotEnergyAngle;
}

// Theme flip — titles/traces bake THEME colors at draw time; replot from cache.
onThemeChange(() => {
    plotEnergyAngle();
    plotMollerXY();
    plotHycalXY();
});
