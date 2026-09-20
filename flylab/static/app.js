/* FlyLab bench UI.
 *
 * No build step: plain ES2020 modules-free script, Plotly + cytoscape from CDN.
 * Chart colour comes from the CSS custom properties in styles.css so light and
 * dark are two selected palettes rather than an automatic flip, and one
 * transmitter palette is shared by every chart and by the circuit viewer.
 */
(function () {
  "use strict";

  // ------------------------------------------------------------------ util
  const $ = (sel, root) => (root || document).querySelector(sel);
  const $$ = (sel, root) => Array.from((root || document).querySelectorAll(sel));

  /* One CDN having a bad minute should not cost the bench a whole panel, so
     each library has a second source that boot() reaches for. These stay
     functions, never captured constants: the fallback lands after this file
     has already been parsed. */
  const hasPlotly = () => typeof window.Plotly !== "undefined";
  const hasCytoscape = () => typeof window.cytoscape !== "undefined";

  const FALLBACKS = {
    Plotly: "https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js",
    cytoscape: "https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.30.2/cytoscape.min.js",
  };

  function loadScript(url) {
    return new Promise((resolve) => {
      const tag = document.createElement("script");
      tag.src = url;
      tag.async = false;
      tag.onload = () => resolve(true);
      tag.onerror = () => resolve(false);
      document.head.appendChild(tag);
    });
  }

  async function ensureLibraries() {
    const jobs = [];
    if (!hasPlotly()) jobs.push(loadScript(FALLBACKS.Plotly));
    if (!hasCytoscape()) jobs.push(loadScript(FALLBACKS.cytoscape));
    if (jobs.length) await Promise.all(jobs);
  }

  function esc(v) {
    return String(v === null || v === undefined ? "" : v).replace(/[&<>"']/g, (c) =>
      ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c])
    );
  }
  function num(v, d) {
    if (v === null || v === undefined || v === "" || Number.isNaN(v)) return "—";
    const n = Number(v);
    if (!Number.isFinite(n)) return "—";
    return n.toFixed(d === undefined ? 3 : d);
  }
  function sci(v, d) {
    const n = Number(v);
    if (!Number.isFinite(n)) return "—";
    if (n === 0) return "0";
    return n.toExponential(d === undefined ? 2 : d);
  }
  function store(key, value) {
    try {
      if (value === undefined) {
        const raw = localStorage.getItem("flylab." + key);
        return raw ? JSON.parse(raw) : null;
      }
      localStorage.setItem("flylab." + key, JSON.stringify(value));
    } catch (err) {
      /* private mode / blocked storage: the bench still works, it just forgets */
    }
    return null;
  }

  // ------------------------------------------------------------------ theme
  function cssVar(name) {
    return getComputedStyle(document.documentElement).getPropertyValue(name).trim() || "#888888";
  }
  function isDark() {
    const stamp = document.documentElement.getAttribute("data-theme");
    if (stamp === "dark") return true;
    if (stamp === "light") return false;
    return window.matchMedia && window.matchMedia("(prefers-color-scheme: dark)").matches;
  }
  function colors() {
    return {
      surface: cssVar("--surface"),
      text1: cssVar("--text-1"),
      text2: cssVar("--text-2"),
      muted: cssVar("--muted"),
      grid: cssVar("--grid"),
      base: cssVar("--baseline"),
      border: cssVar("--border"),
      insect: cssVar("--insect"),
      vertebrate: cssVar("--vertebrate"),
      divLow: cssVar("--div-low"),
      divMid: cssVar("--div-mid"),
      divHigh: cssVar("--div-high"),
      seq: [cssVar("--seq-100"), cssVar("--seq-250"), cssVar("--seq-400"), cssVar("--seq-550"), cssVar("--seq-700")],
      good: cssVar("--good"),
      warning: cssVar("--warning"),
      critical: cssVar("--critical"),
    };
  }
  const NT_ORDER = [
    "acetylcholine",
    "gaba",
    "glutamate",
    "octopamine",
    "dopamine",
    "serotonin",
    "histamine",
    "unclear",
  ];
  function ntColor(nt) {
    const key = NT_ORDER.indexOf(String(nt || "unclear").toLowerCase()) >= 0 ? String(nt).toLowerCase() : "unclear";
    return cssVar("--nt-" + key);
  }
  /* Receptor families get the 7 validated categorical slots in a fixed order,
     so a family keeps its hue when the compound (and the series count) changes. */
  const FAMILY_ORDER = ["nAChR", "GABA_A", "GluCl", "AChE", "Nav", "OctR", "other"];
  function familyColor(family) {
    const slot = FAMILY_ORDER.indexOf(family);
    const vars = [
      "--nt-acetylcholine",
      "--nt-gaba",
      "--nt-glutamate",
      "--nt-octopamine",
      "--nt-dopamine",
      "--nt-serotonin",
      "--nt-histamine",
    ];
    return cssVar(vars[slot < 0 ? vars.length - 1 : Math.min(slot, vars.length - 1)]);
  }
  function sequentialScale() {
    const c = colors().seq;
    return [
      [0, c[0]],
      [0.25, c[1]],
      [0.5, c[2]],
      [0.75, c[3]],
      [1, c[4]],
    ];
  }
  /* Diverging around "no change": blue (down) <-> neutral grey <-> red (up). */
  function divergingColor(delta, scale) {
    const c = colors();
    const t = Math.max(-1, Math.min(1, (delta || 0) / (scale || 1)));
    if (Math.abs(t) < 0.02) return c.divMid;
    return t < 0 ? c.divLow : c.divHigh;
  }

  // ------------------------------------------------------------------ plots
  const PLOT_CFG = { displayModeBar: false, responsive: true };

  function axis(over) {
    const c = colors();
    return Object.assign(
      {
        gridcolor: c.grid,
        zerolinecolor: c.base,
        linecolor: c.base,
        tickcolor: c.base,
        tickfont: { color: c.muted, size: 11 },
        titlefont: { color: c.text2, size: 11.5 },
        automargin: true,
      },
      over || {}
    );
  }
  function layout(over) {
    const c = colors();
    const base = {
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(0,0,0,0)",
      font: { family: 'system-ui, -apple-system, "Segoe UI", Roboto, sans-serif', size: 12, color: c.text2 },
      margin: { t: 14, r: 14, b: 44, l: 58 },
      height: 280,
      bargap: 0.3,
      bargroupgap: 0.12,
      showlegend: false,
      /* Legends sit ABOVE the plot area: a horizontal legend below the x-axis
         wraps at phone width and spills out of its card. */
      legend: {
        orientation: "h",
        x: 0,
        xanchor: "left",
        y: 1.0,
        yanchor: "bottom",
        font: { color: c.text2, size: 11 },
        bgcolor: "rgba(0,0,0,0)",
      },
      hoverlabel: { bgcolor: c.surface, bordercolor: c.base, font: { color: c.text1, size: 11.5 } },
      xaxis: axis(),
      yaxis: axis(),
    };
    const out = Object.assign({}, base, over || {});
    out.xaxis = Object.assign({}, base.xaxis, (over && over.xaxis) || {});
    out.yaxis = Object.assign({}, base.yaxis, (over && over.yaxis) || {});
    out.legend = Object.assign({}, base.legend, (over && over.legend) || {});
    out.margin = Object.assign({}, base.margin, (over && over.margin) || {});
    if (out.showlegend) {
      // room for the legend (two wrapped lines at phone width) plus any title
      out.margin.t = Math.max(out.margin.t || 0, 48);
      out.height = (out.height || 280) + 24;
    }
    return out;
  }
  function draw(id, traces, over) {
    const node = document.getElementById(id);
    if (!node) return;
    if (!hasPlotly()) {
      node.innerHTML = '<div class="empty">Plotly could not be loaded from the CDN; the table view below carries the same values.</div>';
      return;
    }
    /* The container owns the box: Plotly autosizes into an element whose
       height we set here. Passing layout.height instead lets the SVG grow past
       a narrow container and spill into the card below it. */
    const lay = layout(over);
    const height = Math.round(lay.height || 280);
    delete lay.height;
    lay.autosize = true;
    node.style.height = height + "px";
    try {
      window.Plotly.react(node, traces, lay, PLOT_CFG);
    } catch (err) {
      node.innerHTML = '<div class="empty">Chart could not be drawn.</div>';
    }
  }
  function resizePanel(panel) {
    if (!hasPlotly() || !panel) return;
    $$(".plot", panel).forEach((node) => {
      if (node.data) {
        try {
          window.Plotly.Plots.resize(node);
        } catch (err) {
          /* a plot that was never drawn has nothing to resize */
        }
      }
    });
  }

  // ------------------------------------------------------------------ tables
  function table(id, head, rows, opts) {
    const t = document.getElementById(id);
    if (!t) return;
    const o = opts || {};
    t.querySelector("thead").innerHTML =
      "<tr>" + head.map((h) => `<th class="${h.num ? "num" : ""}">${esc(h.label)}</th>`).join("") + "</tr>";
    if (!rows.length) {
      t.querySelector("tbody").innerHTML =
        `<tr><td colspan="${head.length}"><div class="empty">${esc(o.empty || "Nothing to show yet.")}</div></td></tr>`;
      return;
    }
    t.querySelector("tbody").innerHTML = rows
      .map((cells, i) => {
        const attrs = o.rowAttrs ? o.rowAttrs(i) : "";
        return (
          `<tr ${attrs}>` +
          cells.map((cell, j) => `<td class="${head[j] && head[j].num ? "num" : ""}">${cell}</td>`).join("") +
          "</tr>"
        );
      })
      .join("");
  }
  function tile(k, v, u) {
    return `<div class="tile"><div class="k">${esc(k)}</div><div class="v">${v}</div><div class="u">${esc(u || "")}</div></div>`;
  }
  function notices(id, list, kind) {
    const node = document.getElementById(id);
    if (!node) return;
    node.innerHTML = (list || [])
      .map((w) => `<div class="notice ${kind || ""}">${esc(w)}</div>`)
      .join("");
  }

  // ------------------------------------------------------------------ status
  const statusDot = $("#status-dot");
  const statusText = $("#status-text");
  const statusRuntime = $("#status-runtime");

  function status(state, text, runtimeMs) {
    if (statusDot) statusDot.dataset.state = state;
    if (statusText) statusText.textContent = text;
    if (statusRuntime) {
      statusRuntime.textContent =
        runtimeMs === undefined || runtimeMs === null ? "" : `${Math.round(runtimeMs)} ms`;
    }
  }
  function toast(message, kind) {
    const host = $("#toasts");
    if (!host) return;
    const div = document.createElement("div");
    div.className = "toast " + (kind || "");
    div.innerHTML = `<b>${esc(kind === "ok" ? "done" : kind === "info" ? "note" : "error")}</b>${esc(message)}`;
    host.appendChild(div);
    setTimeout(() => div.remove(), kind === "error" || !kind ? 9000 : 5000);
  }

  // ------------------------------------------------------------------ api
  async function api(path, options) {
    const opts = Object.assign({ headers: { Accept: "application/json" } }, options || {});
    let res;
    try {
      res = await fetch(path, opts);
    } catch (err) {
      throw Object.assign(new Error("the bench server is not reachable"), { status: 0 });
    }
    if (res.status === 501) {
      const body = await res.json().catch(() => ({}));
      throw Object.assign(new Error(body.detail || "module not available yet"), { status: 501 });
    }
    if (!res.ok) {
      const body = await res.json().catch(() => ({}));
      const detail = typeof body.detail === "string" ? body.detail : `HTTP ${res.status}`;
      throw Object.assign(new Error(detail), { status: res.status });
    }
    const ctype = res.headers.get("content-type") || "";
    return ctype.indexOf("application/json") >= 0 ? res.json() : res.text();
  }
  const post = (path, body) =>
    api(path, {
      method: "POST",
      headers: { "Content-Type": "application/json", Accept: "application/json" },
      body: JSON.stringify(body),
    });

  // ------------------------------------------------------------------ state
  const state = {
    meta: null,
    census: null,
    notebook: null,
    cache: {},
    loaded: {},
    dirty: {},
    experiment: null,
    runOrder: null,
    blindKey: null,
    revealed: false,
    cy: null,
    activeTab: "scorecard",
  };

  function design() {
    const conc = Number($("#f-conc").value);
    return {
      compound: $("#f-compound").value || null,
      conc_M: Number.isFinite(conc) && conc >= 0 ? conc : 0,
      sugar_hz: Number($("#f-sugar").value) || 0,
      bitter_hz: Number($("#f-bitter").value) || 0,
      engine: $("#f-engine").value,
      graph: $("#f-graph").value,
      seed: Math.max(0, parseInt($("#f-seed").value, 10) || 0),
      genotype: $("#genotype-field").hidden ? null : $("#f-genotype").value || null,
    };
  }

  function setNotebook(nb, label) {
    if (!nb || typeof nb !== "object") return;
    state.notebook = nb;
    renderWarnings();
    renderNotebook();
    pushHistory(nb, label);
  }

  function renderWarnings() {
    const list = (state.notebook && state.notebook.warnings) || [];
    const count = $("#warn-count");
    if (count) count.textContent = String(list.length);
    const body = $("#drawer-body");
    if (body) {
      body.innerHTML = list.length
        ? "<ul>" + list.map((w) => `<li>${esc(w)}</li>`).join("") + "</ul>"
        : '<p class="empty">This run produced no warnings.</p>';
    }
    const nbw = $("#nb-warnings");
    if (nbw) {
      nbw.innerHTML = list.length
        ? list.map((w) => `<div class="notice">${esc(w)}</div>`).join("")
        : '<div class="empty">This run produced no warnings.</div>';
    }
  }

  // ------------------------------------------------------------------ history
  /* Rasters and PSTHs are megabytes; history keeps the notebook minus those two
     so 50 runs fit in the browser's storage quota. The full notebook stays in
     memory for the current run and in any file you export. */
  function forHistory(nb) {
    const copy = JSON.parse(JSON.stringify(nb));
    if (copy.readouts) {
      ["raster", "psth", "top_changed", "top_relays"].forEach((k) => {
        if (copy.readouts[k] !== undefined) copy.readouts[k] = "[trimmed for session history]";
      });
    }
    if (copy.uncertainty && copy.uncertainty.replicates) {
      copy.uncertainty.replicates = "[trimmed for session history]";
    }
    return copy;
  }

  function pushHistory(nb, label) {
    const hist = store("history") || [];
    hist.unshift({
      at: new Date().toISOString(),
      label: label || nb.assay || "run",
      assay: nb.assay || null,
      compound: nb.compound || null,
      conc_M: nb.concentration_M === undefined ? null : nb.concentration_M,
      notebook: forHistory(nb),
    });
    let kept = hist.slice(0, 50);
    while (kept.length) {
      try {
        localStorage.setItem("flylab.history", JSON.stringify(kept));
        break;
      } catch (err) {
        kept = kept.slice(0, Math.floor(kept.length / 2)); // quota: keep the newest half
      }
    }
    renderHistory();
  }
  function renderHistory() {
    const hist = store("history") || [];
    table(
      "tbl-history",
      [{ label: "when" }, { label: "run" }, { label: "compound" }, { label: "conc (M)", num: true }],
      hist.map((h) => [
        esc(h.at.slice(11, 19)),
        esc(h.label),
        esc(h.compound || "vehicle"),
        h.conc_M === null ? "—" : sci(h.conc_M),
      ]),
      { empty: "No runs recorded in this browser yet.", rowAttrs: (i) => `class="clickable" data-hist="${i}"` }
    );
  }

  // ------------------------------------------------------------------ meta
  function receptorMeta(name) {
    const meta = state.meta || {};
    const receptors = meta.receptors || {};
    const pairs = meta.selectivity_pairs || {};
    let family = null;
    Object.keys(pairs).forEach((p) => {
      if (pairs[p].insect === name || pairs[p].vertebrate === name) family = p;
    });
    if (!family) {
      if (/nAChR/i.test(name)) family = "nAChR";
      else if (/RDL|GABA/i.test(name)) family = "GABA_A";
      else if (/GluCl|GlyR/i.test(name)) family = "GluCl";
      else if (/AChE/i.test(name)) family = "AChE";
      else if (/Nav/i.test(name)) family = "Nav";
      else if (/Oct/i.test(name)) family = "OctR";
      else family = "other";
    }
    const organism =
      (receptors[name] && receptors[name].organism) || (name.indexOf("insect") === 0 ? "insect" : "vertebrate");
    return { family: family, organism: organism };
  }

  async function loadMeta() {
    state.meta = await api("/api/meta");
    const meta = state.meta;
    const sel = $("#f-compound");
    const multi = $("#f-exp-compounds");
    sel.innerHTML = meta.compounds
      .map((c) => `<option value="${esc(c.key)}">${esc(c.name)}</option>`)
      .join("");
    multi.innerHTML = meta.compounds
      .map((c) => `<option value="${esc(c.key)}">${esc(c.name)}</option>`)
      .join("");
    const saved = store("design") || {};
    if (saved.compound && meta.compounds.some((c) => c.key === saved.compound)) sel.value = saved.compound;
    else if (meta.compounds.some((c) => c.key === "imidacloprid")) sel.value = "imidacloprid";
    if (saved.graph) $("#f-graph").value = saved.graph;
    if (saved.engine) $("#f-engine").value = saved.engine;
    if (saved.conc_M) setConc(Number(saved.conc_M));

    $("#status-version").textContent = "v" + meta.version;
    $("#status-hash").textContent = "lib " + String(meta.library.sha256 || "").slice(0, 12);
    updateGraphHint();
    updateCompoundBadges();

    try {
      const genos = await api("/api/genotypes");
      const list = Array.isArray(genos) ? genos : genos.genotypes || [];
      if (list.length) {
        $("#f-genotype").innerHTML =
          '<option value="">wild type (library default)</option>' +
          list
            .map((g) => {
              const key = typeof g === "string" ? g : g.id || g.key || g.name;
              const label =
                typeof g === "string"
                  ? g
                  : [g.gene, g.allele].filter(Boolean).join(" ") || g.name || g.id;
              return `<option value="${esc(key)}">${esc(label)}</option>`;
            })
            .join("");
        $("#genotype-field").hidden = false;
      }
    } catch (err) {
      /* genotype module has not landed: the control stays hidden, by design */
    }
  }

  function updateGraphHint() {
    const meta = state.meta;
    if (!meta) return;
    const g = (meta.graphs || {})[$("#f-graph").value];
    const hint = $("#graph-hint");
    if (!g || !g.available) {
      hint.textContent = "graph file not present in this checkout";
      return;
    }
    hint.textContent = `${g.n_nodes} cells, ${g.n_edges} edges, ${g.n_seed_cells} seed cells`;
  }

  function updateCompoundBadges() {
    const meta = state.meta;
    const host = $("#compound-badges");
    if (!meta || !host) return;
    const key = $("#f-compound").value;
    const c = (meta.compounds || []).find((x) => x.key === key);
    if (!c) {
      host.innerHTML = "";
      return;
    }
    const targets = c.insect_targets || [];
    const tier = targets.length ? targets[0].evidence_tier : "class_placeholder";
    host.innerHTML =
      `<span class="badge">${esc(c.class || "unclassified")}</span> ` +
      `<span class="badge tier-${esc(tier)}">${esc(tier.replace(/_/g, " "))}</span>` +
      (targets.length
        ? `<div class="hint">insect targets: ${targets.map((t) => esc(t.receptor.replace("insect_", "")) + " " + sci(t.ec50_M, 1) + " M").join(", ")}</div>`
        : '<div class="hint">no sourced insect target in the library</div>');
  }

  function setConc(value) {
    const v = Number(value);
    if (!Number.isFinite(v) || v <= 0) return;
    $("#f-conc").value = v.toExponential(2).replace("e+", "e");
    $("#f-logconc").value = String(Math.max(-11, Math.min(-3, Math.log10(v))));
    $("#conc-hint").textContent = sci(v) + " M";
  }

  // ==================================================================
  // 1. SCORECARD
  // ==================================================================
  async function runScorecard() {
    const d = design();
    const [occ, nb] = await Promise.all([
      api(`/api/occupancy?compound=${encodeURIComponent(d.compound)}&conc_M=${d.conc_M}`),
      post("/api/assay/subgraph", {
        compound: d.compound,
        conc_M: d.conc_M,
        graph: d.graph,
        genotype: d.genotype,
      }),
    ]);
    state.cache.scorecard = { occ: occ, nb: nb };
    setNotebook(nb, "scorecard · subgraph");
    renderScorecard();
    return nb.readouts && nb.readouts.runtime_ms;
  }

  function renderScorecard() {
    const data = state.cache.scorecard;
    if (!data) return;
    const { occ, nb } = data;
    const c = colors();
    const pairs = occ.selectivity || {};
    const names = Object.keys(pairs);

    const hatch = (flag) => (flag ? "/" : "");
    draw(
      "plot-scorecard",
      [
        {
          type: "bar",
          name: "insect",
          x: names,
          y: names.map((n) => pairs[n].insect_occupancy),
          marker: {
            color: c.insect,
            pattern: { shape: names.map((n) => hatch(pairs[n].placeholder)), fgcolor: c.surface, size: 5, solidity: 0.35 },
            line: { width: 2, color: c.surface },
          },
          hovertemplate: "%{x} insect<br>occupancy %{y:.3f}<extra></extra>",
        },
        {
          type: "bar",
          name: "vertebrate",
          x: names,
          y: names.map((n) => pairs[n].vertebrate_occupancy),
          marker: {
            color: c.vertebrate,
            pattern: { shape: names.map((n) => hatch(pairs[n].placeholder)), fgcolor: c.surface, size: 5, solidity: 0.35 },
            line: { width: 2, color: c.surface },
          },
          hovertemplate: "%{x} vertebrate<br>occupancy %{y:.3f}<extra></extra>",
        },
      ],
      {
        barmode: "group",
        showlegend: true,
        height: 300,
        yaxis: axis({ title: "fractional occupancy", range: [0, 1.02] }),
        xaxis: axis({ title: "" }),
      }
    );
    $("#legend-scorecard").innerHTML =
      `<span class="item" style="color:${c.insect}"><span class="swatch" style="background:${c.insect}"></span>insect target</span>` +
      `<span class="item" style="color:${c.vertebrate}"><span class="swatch" style="background:${c.vertebrate}"></span>vertebrate counterpart</span>` +
      `<span class="item">hatched = class placeholder (no sourced EC50)</span>`;

    table(
      "tbl-selectivity",
      [
        { label: "pair" },
        { label: "insect EC50", num: true },
        { label: "vert EC50", num: true },
        { label: "log10 ratio", num: true },
        { label: "Δ occupancy", num: true },
        { label: "tier" },
      ],
      names.map((n) => {
        const p = pairs[n];
        return [
          esc(n),
          sci(p.insect_ec50_M),
          sci(p.vertebrate_ec50_M),
          num(p.log10_ec50_ratio_vert_over_insect, 2),
          num(p.occupancy_difference, 3),
          `<span class="badge tier-${esc(p.evidence_tier)}">${esc(p.evidence_tier.replace(/_/g, " "))}</span>`,
        ];
      })
    );

    table(
      "tbl-occupancy",
      [
        { label: "receptor" },
        { label: "occupancy", num: true },
        { label: "EC50 (M)", num: true },
        { label: "n", num: true },
        { label: "direction" },
        { label: "tier" },
      ],
      (occ.receptors || []).map((r) => {
        const m = receptorMeta(r.receptor);
        return [
          `<span class="swatch" style="background:${m.organism === "insect" ? c.insect : c.vertebrate}"></span>${esc(r.receptor)}`,
          num(r.occupancy, 3),
          sci(r.ec50_M),
          num(r.n, 1),
          esc(r.direction),
          `<span class="badge tier-${esc(r.evidence_tier)}">${esc(String(r.evidence_tier).replace(/_/g, " "))}</span>`,
        ];
      })
    );

    renderGains(nb.gains || {});
    renderMechanisms(occ.receptors || []);
  }

  function renderGains(gains) {
    const host = $("#meters-gains");
    if (!host) return;
    const keys = ["g_ach", "g_gaba", "g_glu", "g_oct", "g_nav", "ach_tone"];
    const c = colors();
    host.innerHTML = keys
      .map((k) => {
        const v = Number(gains[k]);
        const value = Number.isFinite(v) ? v : 1;
        /* 0 -> 0%, 1 -> 50%, 3 -> 100%: vehicle sits at the track's midpoint. */
        const pos = value <= 1 ? (value / 1) * 50 : 50 + Math.min(1, (value - 1) / 2) * 50;
        const left = Math.min(50, pos);
        const width = Math.max(1.2, Math.abs(pos - 50));
        const fill = Math.abs(value - 1) < 1e-9 ? c.divMid : value < 1 ? c.divLow : c.divHigh;
        return (
          `<div class="meter"><span class="name">${esc(k)}</span>` +
          `<span class="track"><span class="fill" style="left:${left}%;width:${width}%;background:${fill}"></span></span>` +
          `<span class="val">${num(value, 3)}</span></div>`
        );
      })
      .join("");
  }

  function renderMechanisms(rows) {
    const meta = state.meta || {};
    const active = new Set();
    (rows || []).forEach((r) => {
      if (r.direction && r.direction !== "none" && Number(r.occupancy) > 1e-6) {
        active.add(`${r.receptor}:${r.direction}`);
      }
    });
    const all = (meta.mechanisms || []).slice();
    all.sort((a, b) => {
      const ka = active.has(`${a.receptor}:${a.direction}`) ? 0 : 1;
      const kb = active.has(`${b.receptor}:${b.direction}`) ? 0 : 1;
      return ka - kb;
    });
    table(
      "tbl-mechanism",
      [{ label: "" }, { label: "receptor" }, { label: "direction" }, { label: "gain" }, { label: "formula" }],
      all.map((m) => [
        active.has(`${m.receptor}:${m.direction}`)
          ? `<span class="badge tier-literature_order">active</span>`
          : `<span class="badge">—</span>`,
        esc(m.receptor),
        esc(m.direction),
        `<span class="mono">${esc(m.gain)}</span>`,
        `<span class="mono">${esc(m.formula)}</span>`,
      ])
    );
  }

  // ==================================================================
  // 2. CURVES
  // ==================================================================
  async function runCurves() {
    const d = design();
    const curve = await api(`/api/occupancy/curve?compound=${encodeURIComponent(d.compound)}`);
    state.cache.curves = Object.assign({}, state.cache.curves, { curve: curve, conc: d.conc_M });
    renderOccCurve();
    if (!state.cache.curves.ic50) await runIc50();
    else renderIc50();
    return null;
  }

  async function runIc50() {
    const d = design();
    const assay = $("#f-ic50-assay").value;
    const readout = $("#f-ic50-readout").value;
    const body = {
      assay: assay,
      compound: d.compound,
      readout: readout,
      n_boot: 200,
      n_rep: 3,
      seed: d.seed,
      graph: assay === "taste" ? null : d.graph,
      genotype: d.genotype,
    };
    const out = await post("/api/analysis/ic50", body);
    state.cache.curves = Object.assign({}, state.cache.curves, { ic50: out });
    renderIc50();
    return out.runtime_ms;
  }

  function renderOccCurve() {
    const data = state.cache.curves;
    if (!data || !data.curve) return;
    const pts = data.curve.points || [];
    if (!pts.length) return;
    const receptors = Object.keys(pts[0].receptors || {});
    const x = pts.map((p) => Math.log10(p.conc_M));
    const traces = receptors.map((r) => {
      const m = receptorMeta(r);
      return {
        type: "scatter",
        mode: "lines",
        name: r,
        x: x,
        y: pts.map((p) => p.receptors[r]),
        line: {
          color: familyColor(m.family),
          width: 2,
          dash: m.organism === "insect" ? "solid" : "dash",
        },
        hovertemplate: `${esc(r)}<br>10^%{x:.1f} M &middot; occupancy %{y:.3f}<extra></extra>`,
      };
    });
    const c = colors();
    const logc = Number.isFinite(Math.log10(data.conc)) ? Math.log10(data.conc) : -6;
    draw("plot-occcurve", traces, {
      showlegend: true,
      height: 340,
      xaxis: axis({ title: "log10 concentration (M)", range: [-11, -3] }),
      yaxis: axis({ title: "fractional occupancy", range: [0, 1.02] }),
      shapes: [
        {
          type: "line",
          x0: logc,
          x1: logc,
          y0: 0,
          y1: 1.02,
          line: { color: c.muted, width: 2 },
        },
      ],
      annotations: [
        {
          x: logc,
          y: 0.99,
          yref: "paper",
          text: "current dose",
          showarrow: false,
          yanchor: "top",
          xanchor: "left",
          font: { color: c.muted, size: 10.5 },
        },
      ],
      margin: { t: 24, r: 14, b: 44, l: 58 },
    });
  }

  function renderIc50() {
    const out = state.cache.curves && state.cache.curves.ic50;
    if (!out) return;
    const c = colors();
    const pts = out.points || [];
    const x = pts.map((p) => Math.log10(p.conc_M));
    const mean = pts.map((p) => (p.mean === null ? p.value : p.mean));
    const sd = pts.map((p) => p.sd || 0);
    const fit = out.fit || {};
    const ci = (out.ci || {}).ic50;
    const shapes = [];
    const annotations = [];
    if (Number.isFinite(Number(fit.ic50)) && Number(fit.ic50) > 0) {
      const lx = Math.log10(Number(fit.ic50));
      if (ci && Number(ci[0]) > 0 && Number(ci[1]) > 0) {
        shapes.push({
          type: "rect",
          xref: "x",
          yref: "paper",
          x0: Math.log10(Number(ci[0])),
          x1: Math.log10(Number(ci[1])),
          y0: 0,
          y1: 1,
          fillcolor: c.divLow,
          opacity: 0.14,
          line: { width: 0 },
        });
      }
      shapes.push({ type: "line", x0: lx, x1: lx, yref: "paper", y0: 0, y1: 1, line: { color: c.muted, width: 2 } });
      annotations.push({
        x: lx,
        y: 0.99,
        yref: "paper",
        text: "model IC50",
        showarrow: false,
        yanchor: "top",
        xanchor: "left",
        font: { color: c.muted, size: 10.5 },
      });
    }
    draw(
      "plot-ic50",
      [
        {
          type: "scatter",
          mode: "lines+markers",
          name: out.readout,
          x: x,
          y: mean,
          error_y: { type: "data", array: sd, color: c.muted, thickness: 1.5, width: 4 },
          line: { color: c.insect, width: 2 },
          marker: { size: 9, color: c.insect, line: { width: 2, color: c.surface } },
          hovertemplate: `10^%{x:.1f} M<br>${esc(out.readout)} %{y:.3f}<extra></extra>`,
        },
      ],
      {
        height: 300,
        showlegend: false,
        xaxis: axis({ title: "log10 concentration (M)" }),
        yaxis: axis({ title: out.readout + " (model units)" }),
        shapes: shapes,
        annotations: annotations,
        margin: { t: 24, r: 14, b: 44, l: 58 },
      }
    );
    $("#tiles-ic50").innerHTML =
      tile("model IC50", Number.isFinite(Number(fit.ic50)) ? sci(fit.ic50) : "—", "M · model-derived") +
      tile("bootstrap CI", ci ? `${sci(ci[0], 1)} – ${sci(ci[1], 1)}` : "—", "M") +
      tile("Hill slope", num(fit.slope, 2), "") +
      tile("fit r²", num(fit.r2, 3), "") +
      tile("assay", esc(out.assay), out.readout);
    const caption = $("#ic50-caption");
    if (caption) {
      caption.innerHTML =
        "<b>Model-derived, not an animal IC50.</b> It is the dose at which this simulated readout sits halfway between its low- and high-dose plateau, and it follows entirely from the teaching EC50 library and the gain patch rules.";
    }
  }

  // ==================================================================
  // 3. CIRCUIT
  // ==================================================================
  async function runCircuit() {
    const d = design();
    const minw = Math.max(1, Number($("#f-cy-minw").value) || 5);
    const maxe = Math.max(100, Number($("#f-cy-maxe").value) || 1500);
    const q = new URLSearchParams({
      graph: d.graph,
      min_weight: String(minw),
      max_edges: String(maxe),
      conc_M: String(d.conc_M),
    });
    if (d.compound) q.set("compound", d.compound);
    const [viewer, impact] = await Promise.all([
      api("/api/graph?" + q.toString()),
      post("/api/graph/impact", {
        compound: d.compound,
        conc_M: d.conc_M,
        graph: d.graph,
        top_nodes: 20,
        top_edges: 20,
        top_paths: 10,
      }),
    ]);
    state.cache.circuit = { viewer: viewer, impact: impact };
    renderCircuit();
    return viewer.runtime_ms;
  }

  /* The server lays the graph out in layers: one x per hop distance, every cell
     in that hop stacked on one y line. A 1000-cell layer is then 26,000 px tall
     and 0 px wide, which fits the viewport as a vertical sliver. Packing each
     layer into a block keeps the left-to-right hop order and the server's
     within-layer ordering (transmitter, then degree) while giving the viewer a
     shape it can actually draw. */
  function packLayers(nodes) {
    const SPACING = 26;
    const LAYER_GAP = 110;
    const byColumn = new Map();
    nodes.forEach((n) => {
      const key = Math.round(n.x);
      if (!byColumn.has(key)) byColumn.set(key, []);
      byColumn.get(key).push(n);
    });
    const columns = Array.from(byColumn.keys()).sort((a, b) => a - b);
    const pos = {};
    let cursor = 0;
    columns.forEach((key) => {
      const members = byColumn.get(key);
      // a small layer (the seed cells) gets one roomy line so its labels fit
      const small = members.length <= 8;
      const step = small ? 70 : SPACING;
      const wide = small ? 1 : Math.max(1, Math.ceil(Math.sqrt(members.length * 0.55)));
      const tall = Math.ceil(members.length / wide);
      members.forEach((n, i) => {
        pos[n.id] = {
          x: cursor + (i % wide) * step,
          y: (Math.floor(i / wide) - (tall - 1) / 2) * step,
        };
      });
      cursor += (wide - 1) * step + LAYER_GAP;
    });
    return pos;
  }

  function renderCircuit() {
    const data = state.cache.circuit;
    if (!data) return;
    const { viewer, impact } = data;
    const c = colors();
    const mode = $("#f-cy-state").value;

    $("#cy-count").textContent = `${viewer.n_nodes} cells · ${viewer.n_edges} edges`;

    const seedLabel = {};
    Object.keys(viewer.seeds || {}).forEach((t) => {
      (viewer.seeds[t] || []).forEach((b) => {
        seedLabel[b] = t;
      });
    });

    const rates = viewer.nodes.map((n) => (mode === "vehicle" ? n.rate_vehicle : n.rate) || 0);
    const maxRate = Math.max.apply(null, rates.concat([1]));
    const maxDelta = Math.max.apply(
      null,
      viewer.edges.map((e) => Math.abs(e.eff_weight_delta || 0)).concat([1e-6])
    );

    const positions = packLayers(viewer.nodes);
    const cyHost = document.getElementById("cy");
    if (!hasCytoscape()) {
      cyHost.innerHTML =
        '<div class="empty">cytoscape.js could not be loaded from the CDN. The tables below carry the same information.</div>';
    } else {
      const elements = [];
      viewer.nodes.forEach((n, i) => {
        const rate = rates[i];
        elements.push({
          data: {
            id: "n" + n.id,
            label: seedLabel[n.id] ? `${seedLabel[n.id]}` : "",
            size: 8 + 26 * Math.sqrt(Math.max(0, rate) / maxRate),
            color: ntColor(n.nt),
            seed: n.is_seed ? 1 : 0,
            title: `${n.type || n.id} · ${n.nt} · ${rate.toFixed(2)} Hz`,
          },
          position: positions[n.id] || { x: n.x, y: n.y },
        });
      });
      const present = new Set(viewer.nodes.map((n) => "n" + n.id));
      viewer.edges.forEach((e, i) => {
        if (!present.has("n" + e.source) || !present.has("n" + e.target)) return;
        const delta = mode === "vehicle" ? 0 : e.eff_weight_delta || 0;
        elements.push({
          data: {
            id: "e" + i,
            source: "n" + e.source,
            target: "n" + e.target,
            width: 0.5 + 2.5 * Math.min(1, e.weight / 60),
            color: divergingColor(delta, maxDelta),
            opacity: Math.abs(delta) > maxDelta * 0.02 ? 0.85 : 0.35,
          },
        });
      });

      try {
        if (state.cy) {
          state.cy.destroy();
          state.cy = null;
        }
        state.cy = window.cytoscape({
          container: cyHost,
          elements: elements,
          layout: { name: "preset", fit: true, padding: 24 },
          pixelRatio: 1,
          textureOnViewport: true,
          hideEdgesOnViewport: true,
          style: [
            {
              selector: "node",
              style: {
                width: "data(size)",
                height: "data(size)",
                "background-color": "data(color)",
                "border-width": 0,
                label: "data(label)",
                color: c.text1,
                "font-size": 11,
                "text-valign": "top",
                "text-margin-y": -3,
                "text-outline-width": 2.5,
                "text-outline-color": c.surface,
              },
            },
            {
              selector: 'node[seed = 1]',
              style: { "border-width": 2.5, "border-color": c.text1, "z-index": 10 },
            },
            {
              selector: "edge",
              style: {
                width: "data(width)",
                "line-color": "data(color)",
                opacity: "data(opacity)",
                "curve-style": "straight",
                "target-arrow-shape": "none",
              },
            },
            {
              selector: ".highlighted",
              style: { "line-color": c.critical, "border-color": c.critical, width: 4, opacity: 1, "z-index": 20 },
            },
          ],
        });
      } catch (err) {
        cyHost.innerHTML = '<div class="empty">The circuit viewer could not be initialised.</div>';
      }
    }

    $("#legend-cy").innerHTML =
      NT_ORDER.map(
        (nt) => `<span class="item"><span class="swatch" style="background:${ntColor(nt)}"></span>${esc(nt)}</span>`
      ).join("") +
      `<span class="item"><span class="swatch" style="background:${c.divLow}"></span>edge suppressed</span>` +
      `<span class="item"><span class="swatch" style="background:${c.divHigh}"></span>edge enhanced</span>`;

    const nodes = (impact && impact.top_nodes) || [];
    table(
      "tbl-topnodes",
      [
        { label: "bodyId" },
        { label: "type" },
        { label: "transmitter" },
        { label: "vehicle Hz", num: true },
        { label: "treated Hz", num: true },
        { label: "Δ Hz", num: true },
      ],
      nodes.map((n) => [
        `<span class="mono">${esc(n.bodyId)}</span>`,
        esc(n.type || "—"),
        `<span class="swatch" style="background:${ntColor(n.nt)}"></span>${esc(n.nt)}`,
        num(n.rate_vehicle, 2),
        num(n.rate_treated, 2),
        num(n.delta_hz, 2),
      ]),
      { empty: "No node changed at this dose.", rowAttrs: (i) => `class="clickable" data-node="${esc(nodes[i].bodyId)}"` }
    );

    renderPaths(impact);
  }

  function renderPaths(impact) {
    const host = $("#paths-body");
    if (!host) return;
    const paths = impact && impact.grn_to_mn9_paths;
    if (!paths) {
      host.innerHTML =
        '<div class="notice info">Switch the graph to <b>taste_motor</b> to see labellar GRN &rarr; MN9 paths: the named graph has no GRN seeds.</div>';
      return;
    }
    const rows = [];
    ["sweet", "bitter"].forEach((mod) => {
      (paths[mod] || []).slice(0, 8).forEach((p) => rows.push(Object.assign({ modality: mod }, p)));
    });
    state.cache.paths = rows;
    host.innerHTML =
      `<div class="notice">${esc(paths.note || "")}</div>` +
      '<div class="table-wrap"><table id="tbl-paths"><thead></thead><tbody></tbody></table></div>';
    table(
      "tbl-paths",
      [
        { label: "modality" },
        { label: "path" },
        { label: "before", num: true },
        { label: "after", num: true },
        { label: "min syn", num: true },
      ],
      rows.map((p) => [
        `<span class="badge">${esc(p.modality)}</span>`,
        `<span class="mono">${esc((p.path || []).join(" → "))}</span>`,
        num(p.product_before, 2),
        num(p.product_after, 2),
        num(p.min_weight, 0),
      ]),
      { empty: "No GRN → MN9 path within the hop limit.", rowAttrs: (i) => `class="clickable" data-path="${i}"` }
    );
  }

  function highlightPath(index) {
    const rows = state.cache.paths || [];
    const p = rows[index];
    if (!p || !state.cy) return;
    state.cy.elements().removeClass("highlighted");
    (p.path || []).forEach((b) => {
      const n = state.cy.getElementById("n" + b);
      if (n && n.length) n.addClass("highlighted");
    });
    const first = state.cy.getElementById("n" + (p.path || [])[0]);
    if (first && first.length) state.cy.animate({ center: { eles: first }, zoom: 0.8 }, { duration: 260 });
  }

  // ==================================================================
  // 4. SPIKES
  // ==================================================================
  async function runSpikes() {
    const d = design();
    const t_ms = Math.min(3000, Math.max(50, Number($("#f-tms").value) || 500));
    const body = { conc_M: 0, drive_hz: 40, t_ms: t_ms, seed: d.seed, graph: d.graph, genotype: d.genotype };
    const [vehicle, treated] = await Promise.all([
      post("/api/assay/spiking", Object.assign({}, body, { compound: null, conc_M: 0 })),
      post("/api/assay/spiking", Object.assign({}, body, { compound: d.compound, conc_M: d.conc_M })),
    ]);
    state.cache.spikes = { vehicle: vehicle, treated: treated };
    setNotebook(treated, "spiking · treated");
    renderSpikes();
    return (treated.readouts && treated.readouts.runtime_ms) || null;
  }

  function rasterTrace(nb, color) {
    const rows = (nb.readouts && nb.readouts.raster) || [];
    const x = [];
    const y = [];
    const text = [];
    let budget = 6000;
    rows.forEach((row, i) => {
      (row.spikes_ms || []).forEach((t) => {
        if (budget-- <= 0) return;
        x.push(t);
        y.push(i);
        text.push(`${row.type || row.bodyId}${row.is_seed ? " (seed)" : ""}`);
      });
    });
    return {
      trace: {
        type: "scatter",
        mode: "markers",
        x: x,
        y: y,
        text: text,
        marker: { symbol: "line-ns-open", size: 5, color: color, line: { width: 1.2, color: color } },
        hovertemplate: "%{text}<br>%{x:.1f} ms<extra></extra>",
      },
      n: rows.length,
      labels: rows.map((r) => (r.is_seed ? r.type : "")),
    };
  }

  function renderSpikes() {
    const data = state.cache.spikes;
    if (!data) return;
    const c = colors();
    const tMax = (data.treated.readouts && data.treated.readouts.t_ms) || 500;

    [["vehicle", c.muted], ["treated", c.insect]].forEach((pair) => {
      const key = pair[0];
      const built = rasterTrace(data[key], pair[1]);
      const ticks = [];
      const labels = [];
      built.labels.forEach((l, i) => {
        if (l) {
          ticks.push(i);
          labels.push(l);
        }
      });
      draw("plot-raster-" + key, [built.trace], {
        height: 300,
        xaxis: axis({ title: "time (ms)", range: [0, tMax] }),
        yaxis: axis({
          title: "cell",
          tickmode: ticks.length ? "array" : "auto",
          tickvals: ticks,
          ticktext: labels,
          range: [-1, Math.max(built.n, 2)],
        }),
      });
    });

    ["MN9", "DNp01"].forEach((typ) => {
      const id = typ === "MN9" ? "plot-psth-mn9" : "plot-psth-dnp01";
      const traces = [["vehicle", c.muted], ["treated", c.insect]].map((pair) => {
        const psth = (data[pair[0]].readouts && data[pair[0]].readouts.psth) || { cells: {}, bin_ms: 10, n_bins: 0 };
        const cells = psth.cells[typ] || [];
        const bins = psth.n_bins || 0;
        const sum = new Array(bins).fill(0);
        cells.forEach((cell) => (cell.counts || []).forEach((v, i) => (sum[i] += v)));
        return {
          type: "bar",
          name: pair[0],
          x: sum.map((_, i) => i * (psth.bin_ms || 10)),
          y: sum,
          marker: { color: pair[1], line: { width: 2, color: c.surface } },
          hovertemplate: `${pair[0]}<br>%{x} ms &middot; %{y} spikes<extra></extra>`,
        };
      });
      draw(id, traces, {
        barmode: "group",
        showlegend: true,
        height: 260,
        xaxis: axis({ title: "time (ms, 10 ms bins)" }),
        yaxis: axis({ title: "spikes" }),
      });
    });

    const vr = data.vehicle.readouts || {};
    const tr = data.treated.readouts || {};
    table(
      "tbl-spikes",
      [{ label: "readout" }, { label: "vehicle", num: true }, { label: "treated", num: true }, { label: "unit" }],
      [
        ["MN9", num(vr.mn9_hz, 2), num(tr.mn9_hz, 2), "Hz"],
        ["DNp01", num(vr.dnp01_hz, 2), num(tr.dnp01_hz, 2), "Hz"],
        ["network mean", num(vr.mean_hz, 3), num(tr.mean_hz, 3), "Hz"],
        ["active fraction", num(vr.frac_active, 3), num(tr.frac_active, 3), ""],
        ["spikes", String(vr.n_spikes || 0), String(tr.n_spikes || 0), "count"],
        ["runtime", num(vr.runtime_ms, 0), num(tr.runtime_ms, 0), "ms"],
      ]
    );
  }

  // ==================================================================
  // 5. TASTE
  // ==================================================================
  async function runTaste() {
    const d = design();
    const [reduced, mapped] = await Promise.all([
      post("/api/assay/taste", {
        compound: d.compound,
        conc_M: d.conc_M,
        sugar_hz: d.sugar_hz,
        bitter_hz: d.bitter_hz,
        genotype: d.genotype,
      }),
      post("/api/assay/taste-map", {
        compound: d.compound,
        conc_M: d.conc_M,
        sugar_hz: d.sugar_hz,
        bitter_hz: d.bitter_hz,
        engine: d.engine,
        seed: d.seed,
        genotype: d.genotype,
      }),
    ]);
    state.cache.taste = { reduced: reduced, mapped: mapped };
    setNotebook(mapped, "taste · map path");
    renderTaste();
    return (mapped.readouts && mapped.readouts.runtime_ms) || null;
  }

  function renderTaste() {
    const data = state.cache.taste;
    if (!data) return;
    const c = colors();
    const both = ((data.reduced.warnings || []).concat(data.mapped.warnings || [])).filter(
      (w, i, arr) => arr.indexOf(w) === i
    );
    notices("taste-notices", both);

    [["reduced", "plot-taste-reduced", "tiles-taste-reduced"], ["mapped", "plot-taste-map", "tiles-taste-map"]].forEach(
      (spec) => {
        const r = data[spec[0]].readouts || {};
        draw(
          spec[1],
          [
            {
              type: "bar",
              x: ["vehicle + sugar", "drug + sugar", "drug + sugar & bitter"],
              y: [r.mn9_vehicle_sugar_hz, r.mn9_sugar_hz, r.mn9_sugar_bitter_hz],
              marker: { color: [c.muted, c.insect, c.vertebrate], line: { width: 2, color: c.surface } },
              hovertemplate: "%{x}<br>MN9 %{y:.2f} Hz<extra></extra>",
            },
          ],
          {
            height: 250,
            xaxis: axis({ title: "" }),
            yaxis: axis({ title: "MN9 (Hz)" }),
          }
        );
        const veto = r.bitter_veto_ratio;
        document.getElementById(spec[2]).innerHTML =
          tile("MN9 sugar", num(r.mn9_sugar_hz, 2), "Hz") +
          tile("MN9 sugar+bitter", num(r.mn9_sugar_bitter_hz, 2), "Hz") +
          tile("veto ratio", num(veto, 3), veto !== null && veto < 1 ? "bitter suppresses" : "no suppression");
      }
    );

    const a = data.reduced.readouts || {};
    const b = data.mapped.readouts || {};
    table(
      "tbl-taste",
      [{ label: "readout" }, { label: "reduced_taste_v0", num: true }, { label: "taste_motor map", num: true }, { label: "unit" }],
      [
        ["MN9, vehicle + sugar", num(a.mn9_vehicle_sugar_hz, 3), num(b.mn9_vehicle_sugar_hz, 3), "Hz"],
        ["MN9, drug + sugar", num(a.mn9_sugar_hz, 3), num(b.mn9_sugar_hz, 3), "Hz"],
        ["MN9, drug + sugar & bitter", num(a.mn9_sugar_bitter_hz, 3), num(b.mn9_sugar_bitter_hz, 3), "Hz"],
        ["bitter veto ratio", num(a.bitter_veto_ratio, 3), num(b.bitter_veto_ratio, 3), ""],
        ["sweet GRN cells", "— (reduced)", num(b.n_sweet_grn, 0), "cells"],
        ["bitter GRN cells", "— (reduced)", num(b.n_bitter_grn, 0), "cells"],
        ["engine", "rate", esc(b.engine || "rate"), ""],
      ]
    );
  }

  // ==================================================================
  // 6. WHOLE CNS
  // ==================================================================
  async function runWholeCns() {
    const d = design();
    const [nb, census] = await Promise.all([
      post("/api/assay/cns", { compound: d.compound, conc_M: d.conc_M, genotype: d.genotype }),
      state.census ? Promise.resolve(state.census) : api("/api/census"),
    ]);
    state.census = census;
    state.cache.wholecns = { nb: nb, census: census };
    setNotebook(nb, "whole CNS census");
    renderWholeCns();
    return (nb.readouts && nb.readouts.runtime_ms) || null;
  }

  function renderWholeCns() {
    const data = state.cache.wholecns;
    if (!data) return;
    const c = colors();
    const r = data.nb.readouts || {};

    const nts = Object.entries(r.neurotransmitter_counts || {}).sort((a, b) => b[1] - a[1]);
    draw(
      "plot-census-nt",
      [
        {
          type: "bar",
          orientation: "h",
          x: nts.map((e) => e[1]),
          y: nts.map((e) => e[0]),
          marker: { color: nts.map((e) => ntColor(e[0])), line: { width: 2, color: c.surface } },
          hovertemplate: "%{y}<br>%{x} cells<extra></extra>",
        },
      ],
      {
        height: 260,
        margin: { t: 10, r: 14, b: 40, l: 110 },
        xaxis: axis({ title: "traced cells" }),
        yaxis: axis({ title: "", autorange: "reversed" }),
      }
    );

    const scs = Object.entries(r.superclass_counts || {})
      .sort((a, b) => b[1] - a[1])
      .slice(0, 12);
    draw(
      "plot-census-sc",
      [
        {
          type: "bar",
          orientation: "h",
          x: scs.map((e) => e[1]),
          y: scs.map((e) => e[0]),
          marker: { color: c.seq[2], line: { width: 2, color: c.surface } },
          hovertemplate: "%{y}<br>%{x} cells<extra></extra>",
        },
      ],
      {
        height: 300,
        margin: { t: 10, r: 14, b: 40, l: 130 },
        xaxis: axis({ title: "traced cells" }),
        yaxis: axis({ title: "", autorange: "reversed" }),
      }
    );

    const keys = ["cns_excitation_index", "cns_inhibition_index"];
    draw(
      "plot-cns-index",
      [
        {
          type: "bar",
          name: "vehicle",
          x: keys,
          y: keys.map((k) => (r.vehicle || {})[k]),
          marker: { color: c.muted, line: { width: 2, color: c.surface } },
          hovertemplate: "vehicle %{x}<br>%{y:.3f}<extra></extra>",
        },
        {
          type: "bar",
          name: "treated",
          x: keys,
          y: keys.map((k) => (r.treated || {})[k]),
          marker: { color: c.insect, line: { width: 2, color: c.surface } },
          hovertemplate: "treated %{x}<br>%{y:.3f}<extra></extra>",
        },
      ],
      { barmode: "group", showlegend: true, height: 240, yaxis: axis({ title: "index (dimensionless)" }) }
    );
    table(
      "tbl-cns-index",
      [{ label: "index" }, { label: "vehicle", num: true }, { label: "treated", num: true }],
      keys
        .concat(["excitation_inhibition_ratio"])
        .map((k) => [esc(k), num((r.vehicle || {})[k], 3), num((r.treated || {})[k], 3)])
    );

    const gust = (data.census && data.census.gustatory) || {};
    const subs = Object.entries(gust.subclass_counts || {}).sort((a, b) => b[1] - a[1]);
    draw(
      "plot-gustatory",
      [
        {
          type: "bar",
          orientation: "h",
          x: subs.map((e) => e[1]),
          y: subs.map((e) => e[0]),
          marker: { color: c.seq[3], line: { width: 2, color: c.surface } },
          hovertemplate: "%{y}<br>%{x} cells<extra></extra>",
        },
      ],
      {
        height: 260,
        margin: { t: 10, r: 14, b: 40, l: 130 },
        xaxis: axis({ title: `gustatory cells (n = ${gust.n_traced || 0})` }),
        yaxis: axis({ title: "", autorange: "reversed" }),
      }
    );
  }

  // ==================================================================
  // 7. EXPOSURE
  // ==================================================================
  async function runExposure() {
    const d = design();
    const out = await post("/api/exposure", {
      compound: d.compound,
      dose: Number($("#f-dose").value) || 0,
      route: $("#f-route").value,
      t_h: Number($("#f-th").value) || 24,
    });
    state.cache.exposure = out;
    renderExposure();
    return out.runtime_ms;
  }

  function renderExposure() {
    const out = state.cache.exposure;
    if (!out) return;
    const c = colors();
    notices("exposure-notices", out.warnings || []);
    $("#tiles-exposure").innerHTML =
      tile("Cmax", sci(out.cmax), "M") +
      tile("Tmax", num(out.tmax, 2), "h") +
      tile("AUC", sci(out.auc), "M·h") +
      tile("route", esc(out.route), `f=${out.params.f} ka=${out.params.ka_per_h}/h ke=${out.params.ke_per_h}/h`) +
      tile("dose", num(out.dose_nmol, 2), "nmol per fly");

    draw(
      "plot-ct",
      [
        {
          type: "scatter",
          mode: "lines",
          x: out.t_h,
          y: out.conc_M,
          line: { color: c.insect, width: 2 },
          fill: "tozeroy",
          fillcolor: c.insect + "22",
          hovertemplate: "%{x:.2f} h<br>%{y:.3e} M<extra></extra>",
        },
      ],
      { height: 270, xaxis: axis({ title: "time (h)" }), yaxis: axis({ title: "haemolymph C(t) (M)" }) }
    );

    const receptors = Object.keys(out.occupancy || {});
    draw(
      "plot-occt",
      receptors.map((rname) => {
        const m = receptorMeta(rname);
        return {
          type: "scatter",
          mode: "lines",
          name: rname,
          x: out.t_h,
          y: out.occupancy[rname],
          line: { color: familyColor(m.family), width: 2, dash: m.organism === "insect" ? "solid" : "dash" },
          hovertemplate: `${esc(rname)}<br>%{x:.2f} h &middot; %{y:.3f}<extra></extra>`,
        };
      }),
      {
        height: 270,
        showlegend: true,
        xaxis: axis({ title: "time (h)" }),
        yaxis: axis({ title: "fractional occupancy", range: [0, 1.02] }),
      }
    );
  }

  // ==================================================================
  // 8. UNCERTAINTY
  // ==================================================================
  async function runUncertainty() {
    const d = design();
    const n_rep = Math.max(2, Math.min(64, Number($("#f-nrep").value) || 8));
    const [ens, sens] = await Promise.all([
      post("/api/assay/ensemble", {
        assay: "subgraph",
        compound: d.compound,
        conc_M: d.conc_M,
        n_rep: n_rep,
        seed: d.seed,
        graph: d.graph,
        genotype: d.genotype,
      }),
      post("/api/analysis/sensitivity", {
        assay: "subgraph",
        compound: d.compound,
        conc_M: d.conc_M,
        readout: "mean_hz",
        seed: d.seed,
        graph: d.graph,
        genotype: d.genotype,
      }),
    ]);
    state.cache.uncertainty = { ens: ens, sens: sens };
    setNotebook(ens, "ensemble · subgraph");
    renderUncertainty();
    return sens.runtime_ms;
  }

  function renderUncertainty() {
    const data = state.cache.uncertainty;
    if (!data) return;
    const c = colors();
    const u = data.ens.uncertainty || { ci: {}, mean: {}, sd: {} };
    const keys = Object.keys(u.ci || {});
    draw(
      "plot-ensemble",
      [
        {
          type: "scatter",
          mode: "markers",
          x: keys.map((k) => u.mean[k]),
          y: keys,
          error_x: {
            type: "data",
            symmetric: false,
            array: keys.map((k) => u.ci[k][1] - u.mean[k]),
            arrayminus: keys.map((k) => u.mean[k] - u.ci[k][0]),
            color: c.muted,
            thickness: 1.5,
            width: 5,
          },
          marker: { size: 11, color: c.insect, line: { width: 2, color: c.surface } },
          hovertemplate: "%{y}<br>mean %{x:.3f}<extra></extra>",
        },
      ],
      {
        height: 240,
        margin: { t: 10, r: 14, b: 40, l: 100 },
        xaxis: axis({ title: "value (model units)" }),
        yaxis: axis({ title: "" }),
      }
    );
    table(
      "tbl-ensemble",
      [{ label: "readout" }, { label: "mean", num: true }, { label: "sd", num: true }, { label: "2.5%", num: true }, { label: "97.5%", num: true }],
      keys.map((k) => [esc(k), num(u.mean[k], 3), num(u.sd[k], 3), num(u.ci[k][0], 3), num(u.ci[k][1], 3)])
    );

    const rows = (data.sens.rows || []).slice().sort((a, b) => (b.span || 0) - (a.span || 0));
    const base = rows.length ? rows[0].base : 0;
    draw(
      "plot-tornado",
      [
        {
          type: "bar",
          orientation: "h",
          name: "halved",
          y: rows.map((r) => r.param),
          x: rows.map((r) => (r.low || 0) - base),
          marker: { color: c.divLow, line: { width: 2, color: c.surface } },
          hovertemplate: "%{y} halved<br>Δ %{x:.3f}<extra></extra>",
        },
        {
          type: "bar",
          orientation: "h",
          name: "doubled",
          y: rows.map((r) => r.param),
          x: rows.map((r) => (r.high || 0) - base),
          marker: { color: c.divHigh, line: { width: 2, color: c.surface } },
          hovertemplate: "%{y} doubled<br>Δ %{x:.3f}<extra></extra>",
        },
      ],
      {
        barmode: "overlay",
        showlegend: true,
        height: 260,
        margin: { t: 10, r: 14, b: 44, l: 130 },
        xaxis: axis({ title: `change in ${data.sens.readout} vs base (${num(base, 3)})` }),
        yaxis: axis({ title: "" }),
      }
    );
    table(
      "tbl-tornado",
      [{ label: "parameter" }, { label: "halved", num: true }, { label: "base", num: true }, { label: "doubled", num: true }, { label: "span", num: true }],
      rows.map((r) => [esc(r.param), num(r.low, 3), num(r.base, 3), num(r.high, 3), num(r.span, 3)])
    );
  }

  // ==================================================================
  // 9. EXPERIMENT
  // ==================================================================
  function ladder() {
    const start = Number($("#f-exp-start").value);
    const stop = Number($("#f-exp-stop").value);
    const ppd = Math.max(1, Math.min(10, parseInt($("#f-exp-ppd").value, 10) || 1));
    if (!(start > 0) || !(stop > 0)) return [1e-9, 1e-8, 1e-7, 1e-6, 1e-5];
    const a = Math.log10(Math.min(start, stop));
    const b = Math.log10(Math.max(start, stop));
    const steps = Math.max(1, Math.round((b - a) * ppd));
    const out = [];
    for (let i = 0; i <= steps; i++) out.push(Math.pow(10, a + ((b - a) * i) / steps));
    return out.map((v) => Number(v.toPrecision(6)));
  }

  function experimentDesign() {
    const picked = $$("#f-exp-compounds option:checked").map((o) => o.value);
    const readouts = $$("#f-exp-readouts option:checked").map((o) => o.value);
    const d = design();
    return {
      assay: $("#f-exp-assay").value,
      compounds: picked.length ? picked : [d.compound],
      concs_M: ladder(),
      replicates: Math.max(1, Math.min(64, parseInt($("#f-exp-reps").value, 10) || 1)),
      seed: d.seed,
      readouts: readouts.length ? readouts : ["mn9_hz", "mean_hz"],
      include_vehicle: $("#f-exp-vehicle").checked,
      graph: d.graph,
      genotype: d.genotype,
    };
  }

  function updateLadderHint() {
    const l = ladder();
    $("#ladder-hint").textContent = `${l.length} concentrations: ${l.map((v) => sci(v, 1)).join(", ")}`;
  }

  async function runExperiment() {
    const spec = experimentDesign();
    const rows = spec.compounds.length * spec.concs_M.length * spec.replicates;
    if (rows > 2000) throw new Error(`design would make ${rows} rows; the server limit is 2000`);
    store("design", Object.assign({}, design(), { experiment: spec }));
    const out = await post("/api/experiment", spec);
    state.experiment = out;
    state.revealed = false;
    state.runOrder = null;
    buildBlindKey(spec.compounds);
    renderExperiment();
    $("#btn-exp-csv").disabled = false;
    $("#btn-exp-json").disabled = false;
    $("#btn-exp-reveal").disabled = !$("#f-exp-blind").checked;
    return out.runtime_ms;
  }

  function buildBlindKey(compounds) {
    const shuffled = compounds.slice();
    /* deterministic-enough client-side shuffle: the key is exported with the CSV */
    for (let i = shuffled.length - 1; i > 0; i--) {
      const j = Math.floor(Math.random() * (i + 1));
      const tmp = shuffled[i];
      shuffled[i] = shuffled[j];
      shuffled[j] = tmp;
    }
    const key = {};
    shuffled.forEach((c, i) => {
      key[c] = String.fromCharCode(65 + (i % 26)) + (i >= 26 ? String(Math.floor(i / 26)) : "");
    });
    state.blindKey = key;
  }

  function labelFor(compound) {
    if (!compound) return "vehicle";
    if (!$("#f-exp-blind").checked || state.revealed) return compound;
    return (state.blindKey && state.blindKey[compound]) || "?";
  }

  function renderExperiment() {
    const out = state.experiment;
    if (!out) return;
    const c = colors();
    const readouts = (out.design && out.design.readouts) || ["mn9_hz"];
    const primary = readouts[0];
    /* The run order is shuffled once, when the design runs, and then held: a
       re-render (theme switch, reveal) must not reshuffle the bench sheet. */
    let summary = (out.summary || []).slice();
    if ($("#f-exp-randomize").checked) {
      if (!state.runOrder || state.runOrder.length !== summary.length) {
        state.runOrder = summary.map((_, i) => i);
        for (let i = state.runOrder.length - 1; i > 0; i--) {
          const j = Math.floor(Math.random() * (i + 1));
          const tmp = state.runOrder[i];
          state.runOrder[i] = state.runOrder[j];
          state.runOrder[j] = tmp;
        }
      }
      summary = state.runOrder.map((i) => summary[i]);
    }

    table(
      "tbl-experiment",
      [{ label: "compound" }, { label: "conc (M)", num: true }, { label: "n", num: true }].concat(
        readouts.map((r) => ({ label: r + " (mean ± sd)", num: true }))
      ),
      summary.map((row) => {
        return [esc(labelFor(row.compound)), row.compound ? sci(row.conc_M) : "0", String(row.n)].concat(
          readouts.map((r) => {
            const m = row[r + "_mean"];
            const s = row[r + "_sd"];
            return m === null || m === undefined ? "—" : `${num(m, 3)} ± ${num(s, 3)}`;
          })
        );
      }),
      { empty: "Run the design to fill this grid." }
    );

    const compounds = (out.design.compounds || []).slice();
    const concs = (out.design.concs_M || []).slice();
    const z = compounds.map((cmp) =>
      concs.map((cc) => {
        const hit = (out.summary || []).find(
          (s) => s.compound === cmp && Math.abs(Number(s.conc_M) - Number(cc)) < Number(cc) * 1e-6
        );
        return hit ? hit[primary + "_mean"] : null;
      })
    );
    draw(
      "plot-exp-heat",
      [
        {
          type: "heatmap",
          z: z,
          x: concs.map((v) => Math.log10(v)),
          y: compounds.map(labelFor),
          colorscale: sequentialScale(),
          colorbar: { title: { text: primary, side: "right", font: { size: 11, color: c.text2 } }, thickness: 10, outlinewidth: 0 },
          hovertemplate: "%{y}<br>10^%{x:.1f} M<br>" + esc(primary) + " %{z:.3f}<extra></extra>",
          xgap: 2,
          ygap: 2,
        },
      ],
      {
        height: Math.max(220, 60 + compounds.length * 34),
        margin: { t: 10, r: 14, b: 48, l: 110 },
        xaxis: axis({ title: "log10 concentration (M)" }),
        yaxis: axis({ title: "" }),
      }
    );

    $("#exp-status").textContent =
      `${out.n_rows} drug rows + ${(out.vehicle_rows || []).length} vehicle rows` +
      ($("#f-exp-blind").checked && !state.revealed ? " · labels blinded" : "") +
      ($("#f-exp-randomize").checked ? " · order randomised" : "");
  }

  function experimentCsv() {
    const out = state.experiment;
    if (!out) return "";
    let text = out.csv || "";
    if ($("#f-exp-blind").checked) {
      const lines = text.split("\n");
      const key = state.blindKey || {};
      const mapped = lines.map((line, i) => {
        if (i === 0 || !line.trim()) return line;
        const parts = line.split(",");
        if (key[parts[0]]) parts[0] = key[parts[0]];
        return parts.join(",");
      });
      text = mapped.join("\n");
      text +=
        "\n# blinding key (code,compound)\n" +
        Object.keys(key)
          .map((k) => `# ${key[k]},${k}`)
          .join("\n") +
        "\n";
    }
    text += `# FlyLab ${(state.meta && state.meta.version) || ""} simulated data; no row is a measurement from a living fly.\n`;
    return text;
  }

  function download(name, text, mime) {
    const blob = new Blob([text], { type: mime || "text/plain;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = name;
    document.body.appendChild(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 2000);
  }

  // ==================================================================
  // 10. CONTROLS
  // ==================================================================
  async function runControls() {
    const d = design();
    const n = Math.max(1, Math.min(200, Number($("#f-null-n").value) || 5));
    const modes = $$("#f-null-modes option:checked").map((o) => o.value);
    const results = await Promise.allSettled([
      post("/api/analysis/null-panel", {
        compound: d.compound,
        conc_M: d.conc_M,
        n: n,
        seed: d.seed,
        modes: modes.length ? modes : undefined,
        include_taste_map: false,
      }),
      post("/api/analysis/selectivity-landscape", {
        assay: "subgraph",
        graphs: [d.graph],
        compounds: landscapeCompounds(),
      }),
    ]);
    state.cache.controls = { nullPanel: results[0], landscape: results[1] };
    renderControls();
    return null;
  }

  function landscapeCompounds() {
    /* A full library landscape is minutes of circuit time, so the panel asks
       for the current compound plus a small, fixed reference set. */
    const d = design();
    const refs = ["imidacloprid", "fipronil", "nicotine", "abamectin", "permethrin"];
    const known = new Set(((state.meta && state.meta.compounds) || []).map((c) => c.key));
    const out = [d.compound].concat(refs.filter((r) => known.has(r)));
    return out.filter((v, i, arr) => v && arr.indexOf(v) === i).slice(0, 6);
  }

  async function runNullHistogram() {
    const d = design();
    const mode = $("#f-null-hist-mode").value;
    const n = Math.max(2, Math.min(200, Number($("#f-null-n").value) || 5));
    const out = await post("/api/analysis/null", {
      assay: "subgraph",
      compound: d.compound,
      conc_M: d.conc_M,
      mode: mode,
      n: n,
      seed: d.seed,
      graph: d.graph,
    });
    state.cache.nullHist = out;
    renderNullHistogram();
  }

  function renderNullHistogram() {
    const out = state.cache.nullHist;
    const host = document.getElementById("plot-null-hist");
    if (!host) return;
    if (!out) {
      host.innerHTML = '<div class="empty">Pick a mode and draw its null distribution.</div>';
      return;
    }
    const c = colors();
    draw(
      "plot-null-hist",
      [
        {
          type: "histogram",
          name: "null",
          x: out.null_effects || [],
          marker: { color: c.muted, line: { width: 2, color: c.surface } },
          hovertemplate: "null effect %{x:.3f}<extra></extra>",
        },
      ],
      {
        showlegend: true,
        height: 260,
        xaxis: axis({ title: `effect on ${out.readout} (treated - vehicle)` }),
        yaxis: axis({ title: "count" }),
        shapes: [
          {
            type: "line",
            x0: out.real_effect,
            x1: out.real_effect,
            yref: "paper",
            y0: 0,
            y1: 1,
            line: { color: c.divHigh, width: 3 },
          },
        ],
        annotations: [
          {
            x: out.real_effect,
            y: 0.99,
            yref: "paper",
            text: `real effect (z = ${num(out.z, 2)})`,
            showarrow: false,
            yanchor: "top",
            xanchor: "left",
            font: { color: c.divHigh, size: 10.5 },
          },
        ],
        margin: { t: 26, r: 14, b: 44, l: 58 },
      }
    );
  }

  function comingSoon(host, what) {
    host.innerHTML = `<div class="notice info"><b>${esc(what)}</b> comes with the v0.5 analysis layer. The endpoint answered 501 (module not available yet), so the bench is showing this notice instead of failing.</div>`;
  }

  function renderControls() {
    const data = state.cache.controls;
    if (!data) return;
    const c = colors();

    const nullHost = $("#null-body");
    if (data.nullPanel.status !== "fulfilled") {
      comingSoon(nullHost, "The null-model panel");
    } else {
      const panel = data.nullPanel.value;
      const rows = panel.rows || [];
      nullHost.innerHTML =
        '<div id="plot-null" class="plot"></div>' +
        '<div class="table-wrap"><table id="tbl-null"><thead></thead><tbody></tbody></table></div>' +
        '<div class="card-actions" style="margin-top:12px">' +
        '<label for="f-null-hist-mode" class="badge">histogram for</label>' +
        '<select id="f-null-hist-mode" style="width:auto">' +
        rows.map((r) => `<option value="${esc(r.mode)}">${esc(r.mode)}</option>`).join("") +
        "</select>" +
        '<button class="small" id="btn-null-hist" type="button">Draw null distribution</button>' +
        "</div><div id=\"plot-null-hist\" class=\"plot\"></div>" +
        (panel.warnings || []).map((w) => `<div class="notice">${esc(w)}</div>`).join("");

      /* Real effect against each mode's null mean, with the null sd as the bar.
         One axis, one unit: the readout's change from vehicle. */
      draw(
        "plot-null",
        [
          {
            type: "scatter",
            mode: "markers",
            name: "null mean \u00b1 sd",
            x: rows.map((r) => r.null_mean),
            y: rows.map((r) => r.mode),
            error_x: {
              type: "data",
              array: rows.map((r) => r.null_sd || 0),
              color: c.muted,
              thickness: 1.5,
              width: 6,
            },
            marker: { size: 10, color: c.muted, line: { width: 2, color: c.surface } },
            hovertemplate: "%{y}<br>null mean %{x:.3f}<extra></extra>",
          },
          {
            type: "scatter",
            mode: "markers",
            name: "real effect",
            x: rows.map((r) => r.real_effect),
            y: rows.map((r) => r.mode),
            marker: { size: 13, symbol: "diamond", color: c.insect, line: { width: 2, color: c.surface } },
            hovertemplate: "%{y}<br>real effect %{x:.3f}<extra></extra>",
          },
        ],
        {
          showlegend: true,
          height: 260,
          margin: { t: 10, r: 14, b: 48, l: 170 },
          xaxis: axis({ title: `change in ${(rows[0] && rows[0].readout) || "readout"} vs vehicle` }),
          yaxis: axis({ title: "" }),
        }
      );
      table(
        "tbl-null",
        [
          { label: "mode" },
          { label: "assay" },
          { label: "real", num: true },
          { label: "null mean", num: true },
          { label: "null sd", num: true },
          { label: "z", num: true },
          { label: "p (2-sided)", num: true },
        ],
        rows.map((r) => [
          esc(r.mode),
          esc(r.assay || "subgraph"),
          num(r.real_effect, 3),
          num(r.null_mean, 3),
          num(r.null_sd, 3),
          num(r.z, 2),
          num(r.p_two_sided, 4),
        ]),
        { empty: "The panel returned no modes." }
      );
      const btn = $("#btn-null-hist");
      if (btn) {
        btn.addEventListener("click", () => {
          btn.disabled = true;
          status("busy", "drawing null distribution\u2026");
          runNullHistogram()
            .then(() => status("ok", "null distribution ready"))
            .catch((err) => {
              status("error", "null distribution failed");
              toast("null: " + err.message, "error");
            })
            .finally(() => (btn.disabled = false));
        });
      }
      renderNullHistogram();
    }

    const landHost = $("#landscape-body");
    if (data.landscape.status !== "fulfilled") {
      comingSoon(landHost, "The selectivity landscape");
      return;
    }
    const land = data.landscape.value;
    const rows = (land.rows || []).filter((r) => !r.skipped_reason);
    landHost.innerHTML =
      '<div id="plot-landscape" class="plot"></div>' +
      '<div class="table-wrap"><table id="tbl-landscape"><thead></thead><tbody></tbody></table></div>' +
      (land.warnings || []).map((w) => `<div class="notice">${esc(w)}</div>`).join("");
    const graphs = Array.from(new Set(rows.map((r) => r.graph || "named")));
    draw(
      "plot-landscape",
      graphs.map((g, i) => {
        const subset = rows.filter((r) => (r.graph || "named") === g);
        return {
          type: "scatter",
          mode: "markers+text",
          name: g,
          x: subset.map((r) => r.receptor_si_log10),
          y: subset.map((r) => r.circuit_si_log10),
          text: subset.map((r) => r.compound),
          textposition: "top center",
          textfont: { size: 10, color: c.text2 },
          marker: {
            size: 12,
            color: i === 0 ? c.insect : c.vertebrate,
            line: { width: 2, color: c.surface },
          },
          hovertemplate: "%{text}<br>receptor SI %{x:.2f}<br>circuit SI %{y:.2f}<extra></extra>",
        };
      }),
      {
        showlegend: true,
        height: 340,
        xaxis: axis({ title: "receptor selectivity index (log10)" }),
        yaxis: axis({ title: "circuit selectivity index (log10)" }),
      }
    );
    table(
      "tbl-landscape",
      [
        { label: "compound" },
        { label: "graph" },
        { label: "pair" },
        { label: "receptor SI", num: true },
        { label: "circuit SI", num: true },
        { label: "gap", num: true },
        { label: "circuit beats receptor?" },
      ],
      (land.rows || []).map((r) => [
        esc(r.name || r.compound),
        esc(r.graph),
        esc(r.receptor_pair || "\u2014"),
        num(r.receptor_si_log10, 2),
        num(r.circuit_si_log10, 2),
        num(r.si_gap_circuit_minus_receptor, 2),
        r.skipped_reason
          ? `<span class="badge">${esc(r.skipped_reason)}</span>`
          : r.circuit_beats_receptor
          ? '<span class="badge tier-literature_order">yes</span>'
          : '<span class="badge">no</span>',
      ]),
      { empty: "No compound produced a landscape row." }
    );
  }

  // ==================================================================
  // 11. NOTEBOOK
  // ==================================================================
  function renderNotebook() {
    const nb = state.notebook;
    const pre = $("#nb-json");
    if (!nb) {
      if (pre) pre.textContent = "No run yet.";
      return;
    }
    const ordered = Object.assign({ warnings: nb.warnings || [] }, nb);
    if (pre) pre.textContent = JSON.stringify(ordered, null, 2);
    const p = nb.provenance || {};
    table(
      "tbl-provenance",
      [{ label: "field" }, { label: "value" }],
      [
        ["assay", esc(nb.assay)],
        ["created (UTC)", esc(nb.created_utc)],
        ["flylab version", esc(p.flylab_version)],
        ["git sha", `<span class="mono">${esc(p.git_sha || "—")}</span>`],
        ["library sha256", `<span class="mono">${esc(String(p.library_sha256 || "").slice(0, 24))}</span>`],
        ["map", esc((p.map && p.map.name) || (nb.map && nb.map.name))],
        ["map citation", esc((p.map && p.map.citation) || (nb.map && nb.map.citation))],
        ["rng seed", esc(p.rng_seed)],
        ["platform", esc(p.platform)],
        ["live_lab", nb.live_lab ? `${esc(nb.live_lab.assay_name)} (${nb.live_lab.n_rows} rows)` : "null (nothing imported)"],
      ]
    );
  }

  async function importLiveLab() {
    if (!state.notebook) throw new Error("run an assay first: live_lab attaches to a notebook");
    const nb = await post("/api/notebook/live-lab", {
      notebook: state.notebook,
      csv_text: $("#f-ll-csv").value,
      assay_name: $("#f-ll-name").value || "live_lab",
      source: $("#f-ll-source").value || "user import",
    });
    setNotebook(nb, "live_lab import");
    toast(`attached ${nb.live_lab.n_rows} measured rows`, "ok");
  }

  // ==================================================================
  // tabs
  // ==================================================================
  const TABS = {
    scorecard: { run: runScorecard, render: renderScorecard },
    curves: { run: runCurves, render: () => { renderOccCurve(); renderIc50(); } },
    circuit: { run: runCircuit, render: renderCircuit },
    spikes: { run: runSpikes, render: renderSpikes },
    taste: { run: runTaste, render: renderTaste },
    wholecns: { run: runWholeCns, render: renderWholeCns },
    exposure: { run: runExposure, render: renderExposure },
    uncertainty: { run: runUncertainty, render: renderUncertainty },
    experiment: { run: runExperiment, render: renderExperiment },
    controls: { run: runControls, render: renderControls },
    notebook: { run: async () => { renderNotebook(); renderHistory(); return null; }, render: renderNotebook },
  };

  /* Tabs that cost minutes of circuit time never start on their own. */
  const MANUAL_TABS = ["controls"];

  let running = false;
  async function runTab(name) {
    const tab = TABS[name];
    if (!tab || running) return;
    running = true;
    const buttons = $$("button[data-run], #btn-run");
    buttons.forEach((b) => (b.disabled = true));
    status("busy", `running ${name}…`);
    const t0 = performance.now();
    try {
      const serverMs = await tab.run();
      state.loaded[name] = true;
      status("ok", `${name} ready`, serverMs || performance.now() - t0);
    } catch (err) {
      status("error", `${name} failed`);
      if (err && err.status === 501) {
        toast(`${name}: ${err.message}`, "info");
      } else {
        toast(`${name}: ${(err && err.message) || "unknown error"}`, "error");
      }
    } finally {
      running = false;
      buttons.forEach((b) => (b.disabled = false));
      $("#btn-exp-csv").disabled = !state.experiment;
      $("#btn-exp-json").disabled = !state.experiment;
    }
  }

  function activate(name) {
    state.activeTab = name;
    $$('[role="tab"]').forEach((b) => b.setAttribute("aria-selected", String(b.id === "tab-" + name)));
    $$(".panel").forEach((p) => p.classList.toggle("active", p.id === "panel-" + name));
    const panel = document.getElementById("panel-" + name);
    if (state.loaded[name] && TABS[name]) {
      try {
        TABS[name].render();
        delete state.dirty[name];
      } catch (err) {
        /* a stale cache should never break tab switching */
      }
      resizePanel(panel);
      scheduleResize();
    } else if (TABS[name] && MANUAL_TABS.indexOf(name) < 0) {
      runTab(name);
    } else {
      resizePanel(panel);
    }
  }

  /* Only the visible panel is ever drawn: Plotly cannot measure a display:none
     container, and a chart laid out at zero width does not recover on its own.
     Hidden panels are marked dirty and redrawn when they are next shown. */
  function rerenderAll() {
    Object.keys(state.loaded).forEach((name) => {
      if (state.loaded[name] && name !== state.activeTab) state.dirty[name] = true;
    });
    if (state.loaded[state.activeTab] && TABS[state.activeTab]) {
      try {
        TABS[state.activeTab].render();
      } catch (err) {
        /* a theme change must never throw */
      }
    }
    scheduleResize();
  }

  function scheduleResize() {
    const panel = document.getElementById("panel-" + state.activeTab);
    requestAnimationFrame(() => resizePanel(panel));
  }

  function invalidate() {
    state.cache = {};
    state.loaded = {};
    state.dirty = {};
  }

  // ==================================================================
  // wiring
  // ==================================================================
  function wire() {
    $$('[role="tab"]').forEach((b) =>
      b.addEventListener("click", () => activate(b.id.replace("tab-", "")))
    );
    $$("button[data-run]").forEach((b) =>
      b.addEventListener("click", () => runTab(b.getAttribute("data-run")))
    );
    $("#btn-run").addEventListener("click", () => runTab(state.activeTab));

    $("#f-logconc").addEventListener("input", (e) => {
      setConc(Math.pow(10, Number(e.target.value)));
      invalidate();
    });
    $("#f-conc").addEventListener("change", (e) => {
      const v = Number(e.target.value);
      if (Number.isFinite(v) && v > 0) setConc(v);
      else toast("concentration must be a positive number in molar", "error");
      invalidate();
    });
    $("#f-compound").addEventListener("change", () => {
      updateCompoundBadges();
      invalidate();
      store("design", design());
    });
    $("#f-graph").addEventListener("change", () => {
      updateGraphHint();
      $("#f-cy-minw").value = $("#f-graph").value === "taste_motor" ? 12 : 5;
      invalidate();
      store("design", design());
    });
    ["#f-engine", "#f-seed", "#f-sugar", "#f-bitter", "#f-genotype"].forEach((sel) =>
      $(sel).addEventListener("change", () => {
        invalidate();
        store("design", design());
      })
    );

    $("#btn-ic50").addEventListener("click", () =>
      runIc50().catch((err) => toast("IC50: " + err.message, "error"))
    );
    $("#f-cy-state").addEventListener("change", renderCircuit);
    $("#btn-cy-fit").addEventListener("click", () => {
      if (state.cy) state.cy.fit(undefined, 24);
    });

    document.addEventListener("click", (ev) => {
      const pathRow = ev.target.closest && ev.target.closest("[data-path]");
      if (pathRow) {
        $$("#tbl-paths tr").forEach((tr) => tr.classList.remove("selected"));
        pathRow.classList.add("selected");
        highlightPath(Number(pathRow.getAttribute("data-path")));
        return;
      }
      const nodeRow = ev.target.closest && ev.target.closest("[data-node]");
      if (nodeRow && state.cy) {
        const n = state.cy.getElementById("n" + nodeRow.getAttribute("data-node"));
        if (n && n.length) {
          state.cy.elements().removeClass("highlighted");
          n.addClass("highlighted");
          state.cy.animate({ center: { eles: n }, zoom: 1.1 }, { duration: 260 });
        }
        return;
      }
      const histRow = ev.target.closest && ev.target.closest("[data-hist]");
      if (histRow) {
        const hist = store("history") || [];
        const entry = hist[Number(histRow.getAttribute("data-hist"))];
        if (entry) {
          state.notebook = entry.notebook;
          renderWarnings();
          renderNotebook();
          toast("loaded " + entry.label + " from history", "info");
        }
      }
    });

    ["#f-exp-start", "#f-exp-stop", "#f-exp-ppd"].forEach((sel) =>
      $(sel).addEventListener("input", updateLadderHint)
    );
    $("#f-exp-blind").addEventListener("change", () => {
      state.revealed = false;
      $("#btn-exp-reveal").disabled = !$("#f-exp-blind").checked || !state.experiment;
      renderExperiment();
    });
    $("#f-exp-randomize").addEventListener("change", () => {
      state.runOrder = null;
      renderExperiment();
    });
    $("#btn-exp-reveal").addEventListener("click", () => {
      state.revealed = !state.revealed;
      $("#btn-exp-reveal").textContent = state.revealed ? "Hide codes" : "Reveal codes";
      renderExperiment();
    });
    $("#btn-exp-csv").addEventListener("click", () =>
      download("flylab_experiment.csv", experimentCsv(), "text/csv;charset=utf-8")
    );
    $("#btn-exp-json").addEventListener("click", () =>
      download("flylab_experiment.json", JSON.stringify(state.experiment, null, 2), "application/json")
    );

    $("#btn-nb-export").addEventListener("click", () => {
      if (!state.notebook) return toast("nothing to export yet: run an assay first", "info");
      download("flylab_notebook.json", JSON.stringify(state.notebook, null, 2), "application/json");
    });
    $("#btn-ll-import").addEventListener("click", () =>
      importLiveLab().catch((err) => toast("live_lab import: " + err.message, "error"))
    );
    $("#btn-history-clear").addEventListener("click", () => {
      store("history", []);
      renderHistory();
    });

    const drawer = $("#drawer");
    $("#btn-warnings").addEventListener("click", () => {
      const open = drawer.dataset.open === "true";
      drawer.dataset.open = String(!open);
      drawer.setAttribute("aria-hidden", String(open));
      $("#btn-warnings").setAttribute("aria-expanded", String(!open));
    });
    $("#btn-drawer-close").addEventListener("click", () => {
      drawer.dataset.open = "false";
      drawer.setAttribute("aria-hidden", "true");
      $("#btn-warnings").setAttribute("aria-expanded", "false");
    });

    $("#btn-theme").addEventListener("click", () => {
      const next = isDark() ? "light" : "dark";
      document.documentElement.setAttribute("data-theme", next);
      store("theme", next);
      rerenderAll();
    });

    document.addEventListener("keydown", (ev) => {
      if (ev.key !== "r" && ev.key !== "R") return;
      if (ev.metaKey || ev.ctrlKey || ev.altKey) return;
      const t = ev.target;
      if (t && (t.tagName === "INPUT" || t.tagName === "TEXTAREA" || t.tagName === "SELECT")) return;
      ev.preventDefault();
      runTab(state.activeTab);
    });

    let resizeTimer = null;
    window.addEventListener("resize", () => {
      if (resizeTimer) clearTimeout(resizeTimer);
      // re-draw rather than only resize: a chart laid out at another width
      // keeps that width's tick density and legend wrapping otherwise
      resizeTimer = setTimeout(() => {
        if (state.loaded[state.activeTab] && TABS[state.activeTab]) {
          try {
            TABS[state.activeTab].render();
          } catch (err) {
            /* a resize must never throw */
          }
        }
        scheduleResize();
      }, 180);
    });
  }

  // ==================================================================
  // boot
  // ==================================================================
  async function boot() {
    const theme = store("theme");
    if (theme === "light" || theme === "dark") document.documentElement.setAttribute("data-theme", theme);

    wire();
    renderHistory();
    updateLadderHint();

    await ensureLibraries();
    if (!hasPlotly()) toast("Plotly did not load from either CDN; charts fall back to tables", "info");
    if (!hasCytoscape()) toast("cytoscape.js did not load from either CDN; the circuit view falls back to tables", "info");

    try {
      await loadMeta();
    } catch (err) {
      status("error", "could not reach the bench server");
      toast("meta: " + err.message, "error");
      return;
    }

    const saved = (store("design") || {}).experiment;
    if (saved) {
      if (saved.assay) $("#f-exp-assay").value = saved.assay;
      if (saved.replicates) $("#f-exp-reps").value = saved.replicates;
      $$("#f-exp-compounds option").forEach((o) => {
        o.selected = (saved.compounds || []).indexOf(o.value) >= 0;
      });
      updateLadderHint();
    } else {
      $$("#f-exp-compounds option").forEach((o) => {
        o.selected = o.value === $("#f-compound").value || o.value === "nicotine";
      });
    }
    $("#f-exp-reps").value = $("#f-exp-reps").value || 2;

    status("ok", "ready");
    activate("scorecard");
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", boot);
  else boot();
})();
