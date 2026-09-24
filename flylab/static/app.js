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
    // schema v3: a missing value is "not modelled", never 0 (see evidence.py)
    if (v === null || v === undefined || v === "") return NOT_MODELLED;
    const n = Number(v);
    if (!Number.isFinite(n)) return "—";
    if (n === 0) return "0";
    return n.toExponential(d === undefined ? 2 : d);
  }
  // Placeholder rows carry no number at all: show that, do not print 0.000.
  const NOT_MODELLED = '<span class="muted">not modelled</span>';
  function engagement(v, d) {
    if (v === null || v === undefined || v === "") return NOT_MODELLED;
    return num(v, d === undefined ? 3 : d);
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
  /* Two transports, one contract. The served bench talks HTTP to
     flylab/server.py; the static bench (GitHub Pages) talks to
     window.flylabCall, which is flylab.browser.bridge.call running in
     Pyodide. Same routes, same payloads, same JSON, so nothing below this
     function knows or cares which one is in use. */
  const isStatic = () => window.FLYLAB_STATIC === true;

  /* Both transports report failures as {"error": {status, detail}} or as an
     HTTP status; this turns either into the Error the panels already expect. */
  function apiError(status, detail) {
    return Object.assign(new Error(detail || `HTTP ${status}`), { status: status });
  }
  function unwrap(body) {
    if (body && typeof body === "object" && body.error && typeof body.error === "object") {
      throw apiError(body.error.status || 500, body.error.detail);
    }
    return body;
  }

  async function staticApi(path, opts) {
    let payload = null;
    if (opts && opts.body) {
      try {
        payload = JSON.parse(opts.body);
      } catch (err) {
        payload = null;
      }
    }
    let body;
    try {
      body = await window.flylabCall(path, payload);
    } catch (err) {
      throw apiError(0, "the in-browser bench is not running");
    }
    return unwrap(body);
  }

  async function api(path, options) {
    const opts = Object.assign({ headers: { Accept: "application/json" } }, options || {});
    if (isStatic()) return staticApi(path, opts);
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
    return ctype.indexOf("application/json") >= 0 ? unwrap(await res.json()) : res.text();
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
    prov: {},
    activeTab: "dashboard",
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

  /* Keep the primary bench controls and the saved experiment design together.
     Editing one setting must not silently discard the experiment specification. */
  function persistDesign(extra) {
    const saved = store("design") || {};
    store("design", Object.assign({}, saved, design(), extra || {}));
  }

  const PRESET_FORMAT = "flylab.design-preset";
  const PRESET_VERSION = 1;

  function optionHas(select, value) {
    return Array.from(select.options).some((option) => option.value === value);
  }

  function selectedOptions(id, fallback) {
    const values = $$("#" + id + " option:checked").map((option) => option.value);
    return values.length ? values : fallback.slice();
  }

  function experimentSettings() {
    const d = design();
    return {
      assay: $("#f-exp-assay").value,
      compounds: selectedOptions("f-exp-compounds", [d.compound].filter(Boolean)),
      ladder_start_M: Number($("#f-exp-start").value),
      ladder_stop_M: Number($("#f-exp-stop").value),
      points_per_decade: Number($("#f-exp-ppd").value),
      replicates: Number($("#f-exp-reps").value),
      readouts: selectedOptions("f-exp-readouts", ["mn9_hz", "mean_hz"]),
      randomize: $("#f-exp-randomize").checked,
      blind: $("#f-exp-blind").checked,
      include_vehicle: $("#f-exp-vehicle").checked,
    };
  }

  function validateExperimentSettings(value) {
    if (!value || typeof value !== "object" || Array.isArray(value)) throw new Error("experiment settings must be an object");
    const assaySelect = $("#f-exp-assay");
    const compoundsSelect = $("#f-exp-compounds");
    const readoutsSelect = $("#f-exp-readouts");
    if (typeof value.assay !== "string" || !optionHas(assaySelect, value.assay)) {
      throw new Error("experiment assay is not available in this FlyLab build");
    }
    const validList = (items, select, label) => {
      if (!Array.isArray(items) || items.length === 0 || items.some((item) => typeof item !== "string" || !optionHas(select, item))) {
        throw new Error("experiment " + label + " contains an unavailable value");
      }
      if (new Set(items).size !== items.length) throw new Error("experiment " + label + " must not contain duplicates");
      return items.slice();
    };
    const compounds = validList(value.compounds, compoundsSelect, "compounds");
    const readouts = validList(value.readouts, readoutsSelect, "readouts");
    const positive = (v) => typeof v === "number" && Number.isFinite(v) && v > 0;
    if (!positive(value.ladder_start_M) || !positive(value.ladder_stop_M)) {
      throw new Error("concentration ladder start and stop must be positive finite numbers");
    }
    if (!Number.isInteger(value.points_per_decade) || value.points_per_decade < 1 || value.points_per_decade > 10) {
      throw new Error("points per decade must be an integer from 1 to 10");
    }
    if (!Number.isInteger(value.replicates) || value.replicates < 1 || value.replicates > 64) {
      throw new Error("replicates must be an integer from 1 to 64");
    }
    ["randomize", "blind", "include_vehicle"].forEach((key) => {
      if (typeof value[key] !== "boolean") throw new Error("experiment " + key + " must be true or false");
    });
    const decades = Math.abs(Math.log10(value.ladder_stop_M) - Math.log10(value.ladder_start_M));
    const concentrationCount = Math.max(1, Math.round(decades * value.points_per_decade)) + 1;
    if (compounds.length * concentrationCount * value.replicates > 2000) {
      throw new Error("experiment preset exceeds the 2,000-row design limit");
    }
    return {
      assay: value.assay,
      compounds: compounds,
      ladder_start_M: value.ladder_start_M,
      ladder_stop_M: value.ladder_stop_M,
      points_per_decade: value.points_per_decade,
      replicates: value.replicates,
      readouts: readouts,
      randomize: value.randomize,
      blind: value.blind,
      include_vehicle: value.include_vehicle,
    };
  }

  function validatePreset(value) {
    if (!value || typeof value !== "object" || Array.isArray(value)) throw new Error("preset file must contain a JSON object");
    if (value.format !== PRESET_FORMAT) throw new Error("file is not a FlyLab design preset");
    if (value.version !== PRESET_VERSION) throw new Error("unsupported preset version; expected version " + PRESET_VERSION);
    const raw = value.design;
    if (!raw || typeof raw !== "object" || Array.isArray(raw)) throw new Error("preset is missing its main design settings");
    if (typeof raw.compound !== "string" || !optionHas($("#f-compound"), raw.compound)) {
      throw new Error("preset compound is not available in this FlyLab build");
    }
    if (typeof raw.conc_M !== "number" || !Number.isFinite(raw.conc_M) || raw.conc_M <= 0 || raw.conc_M > 1) {
      throw new Error("preset concentration must be greater than 0 and at most 1 M");
    }
    if (!["named", "taste_motor"].includes(raw.graph) || !optionHas($("#f-graph"), raw.graph)) {
      throw new Error("preset graph is not available in this FlyLab build");
    }
    if (!["rate", "lif"].includes(raw.engine) || !optionHas($("#f-engine"), raw.engine)) {
      throw new Error("preset engine is not available in this FlyLab build");
    }
    ["sugar_hz", "bitter_hz"].forEach((key) => {
      if (typeof raw[key] !== "number" || !Number.isFinite(raw[key]) || raw[key] < 0 || raw[key] > 1000) {
        throw new Error(key.replace("_hz", " drive") + " must be between 0 and 1,000 Hz");
      }
    });
    if (!Number.isSafeInteger(raw.seed) || raw.seed < 0 || raw.seed > 2147483647) {
      throw new Error("RNG seed must be an integer from 0 to 2,147,483,647");
    }
    if (raw.genotype !== null && typeof raw.genotype !== "string") throw new Error("genotype must be a string or null");
    if (raw.genotype !== null && (!$("#genotype-field") || $("#genotype-field").hidden || !optionHas($("#f-genotype"), raw.genotype))) {
      throw new Error("preset genotype is not available in this FlyLab build");
    }
    return {
      design: {
        compound: raw.compound,
        conc_M: raw.conc_M,
        sugar_hz: raw.sugar_hz,
        bitter_hz: raw.bitter_hz,
        engine: raw.engine,
        graph: raw.graph,
        seed: raw.seed,
        genotype: raw.genotype,
      },
      experiment: validateExperimentSettings(value.experiment),
    };
  }

  function applyExperimentSettings(settings) {
    $("#f-exp-assay").value = settings.assay;
    $("#f-exp-start").value = String(settings.ladder_start_M);
    $("#f-exp-stop").value = String(settings.ladder_stop_M);
    $("#f-exp-ppd").value = String(settings.points_per_decade);
    $("#f-exp-reps").value = String(settings.replicates);
    $$("#f-exp-compounds option").forEach((option) => {
      option.selected = settings.compounds.includes(option.value);
    });
    $$("#f-exp-readouts option").forEach((option) => {
      option.selected = settings.readouts.includes(option.value);
    });
    $("#f-exp-randomize").checked = settings.randomize;
    $("#f-exp-blind").checked = settings.blind;
    $("#f-exp-vehicle").checked = settings.include_vehicle;
    updateLadderHint();
  }

  function presetForExport() {
    const preset = { format: PRESET_FORMAT, version: PRESET_VERSION, design: design(), experiment: experimentSettings() };
    validatePreset(preset);
    return preset;
  }

  function presetStatus(message) {
    const node = $("#preset-status");
    if (node) node.textContent = message;
  }

  function applyPreset(preset) {
    const validated = validatePreset(preset);
    const main = validated.design;
    $("#f-compound").value = main.compound;
    setConc(main.conc_M);
    $("#f-sugar").value = String(main.sugar_hz);
    $("#f-bitter").value = String(main.bitter_hz);
    $("#f-engine").value = main.engine;
    $("#f-graph").value = main.graph;
    $("#f-seed").value = String(main.seed);
    if (!$("#genotype-field").hidden) $("#f-genotype").value = main.genotype || "";
    updateCompoundBadges();
    updateGraphHint();
    const minw = $("#f-cy-minw");
    if (minw) minw.value = main.graph === "taste_motor" ? 12 : 5;
    applyExperimentSettings(validated.experiment);

    state.experiment = null;
    state.runOrder = null;
    state.blindKey = null;
    state.revealed = false;
    $("#btn-exp-csv").disabled = true;
    $("#btn-exp-json").disabled = true;
    $("#btn-exp-reveal").disabled = true;
    $("#btn-exp-reveal").textContent = "Reveal codes";
    invalidate();
    renderExperiment();
    persistDesign({ experiment: experimentDesign(), experimentSettings: experimentSettings() });
  }

  function restoreExperimentSettings(saved) {
    if (saved.experimentSettings) {
      try {
        applyExperimentSettings(validateExperimentSettings(saved.experimentSettings));
        return;
      } catch (err) {
        /* Invalid or stale browser storage falls back to the older saved design below. */
      }
    }
    const previous = saved.experiment;
    if (previous) {
      if (typeof previous.assay === "string" && optionHas($("#f-exp-assay"), previous.assay)) {
        $("#f-exp-assay").value = previous.assay;
      }
      if (Number.isInteger(previous.replicates) && previous.replicates >= 1 && previous.replicates <= 64) {
        $("#f-exp-reps").value = String(previous.replicates);
      }
      $$("#f-exp-compounds option").forEach((option) => {
        option.selected = Array.isArray(previous.compounds) && previous.compounds.includes(option.value);
      });
      const concentrations = Array.isArray(previous.concs_M) ? previous.concs_M.filter((value) => Number.isFinite(value) && value > 0) : [];
      if (concentrations.length > 1) {
        $("#f-exp-start").value = String(Math.min.apply(null, concentrations));
        $("#f-exp-stop").value = String(Math.max.apply(null, concentrations));
      }
      if (Array.isArray(previous.readouts)) {
        $$("#f-exp-readouts option").forEach((option) => {
          option.selected = previous.readouts.includes(option.value);
        });
      }
      if (typeof previous.include_vehicle === "boolean") $("#f-exp-vehicle").checked = previous.include_vehicle;
      updateLadderHint();
    } else {
      $$("#f-exp-compounds option").forEach((option) => {
        option.selected = option.value === $("#f-compound").value || option.value === "nicotine";
      });
    }
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
    const options = meta.compounds
      .map((c) => `<option value="${esc(c.key)}">${esc(c.name)}</option>`)
      .join("");
    ["#f-cmp-compounds", "#f-mix-a", "#f-mix-b", "#f-geno-other"].forEach((id) => {
      const node = $(id);
      if (node) node.innerHTML = options;
    });
    const pick = (id, key) => {
      const node = $(id);
      if (node && meta.compounds.some((c) => c.key === key)) node.value = key;
    };
    pick("#f-mix-a", "imidacloprid");
    pick("#f-mix-b", "fipronil");
    pick("#f-geno-other", "ddt");
    $$("#f-cmp-compounds option").forEach((o) => {
      o.selected = ["imidacloprid", "fipronil", "deltamethrin"].indexOf(o.value) >= 0;
    });

    const saved = store("design") || {};
    if (saved.compound && meta.compounds.some((c) => c.key === saved.compound)) sel.value = saved.compound;
    else if (meta.compounds.some((c) => c.key === "imidacloprid")) sel.value = "imidacloprid";
    if (["named", "taste_motor"].includes(saved.graph)) $("#f-graph").value = saved.graph;
    if (["rate", "lif"].includes(saved.engine)) $("#f-engine").value = saved.engine;
    if (Number.isFinite(Number(saved.sugar_hz)) && Number(saved.sugar_hz) >= 0) {
      $("#f-sugar").value = Math.min(1000, Number(saved.sugar_hz));
    }
    if (Number.isFinite(Number(saved.bitter_hz)) && Number(saved.bitter_hz) >= 0) {
      $("#f-bitter").value = Math.min(1000, Number(saved.bitter_hz));
    }
    if (Number.isSafeInteger(Number(saved.seed)) && Number(saved.seed) >= 0) {
      $("#f-seed").value = Number(saved.seed);
    }
    if (Number.isFinite(Number(saved.conc_M)) && Number(saved.conc_M) > 0) setConc(Number(saved.conc_M));

    // re-apply the current concentration so the dashboard rail and its ladder
    // buttons show the right value on the first paint, saved design or not
    setConc(Number($("#f-conc").value) || 1e-6);

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
        if (saved.genotype && Array.from($("#f-genotype").options).some((o) => o.value === saved.genotype)) {
          $("#f-genotype").value = saved.genotype;
        }
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

  /* One concentration, three controls: the sidebar slider, the sidebar box and
     the dashboard dose rail all read and write the same value. */
  function setConc(value) {
    const v = Number(value);
    if (!Number.isFinite(v) || v <= 0) return;
    const log = String(Math.max(-11, Math.min(-3, Math.log10(v))));
    $("#f-conc").value = v.toExponential(2).replace("e+", "e");
    $("#f-logconc").value = log;
    $("#conc-hint").textContent = sci(v) + " M";
    const rail = $("#f-dash-logconc");
    if (rail) rail.value = log;
    const read = $("#dash-conc-read");
    if (read) read.textContent = sci(v) + " M";
    $$("#dash-ladder button").forEach((b) => {
      const step = Number(b.getAttribute("data-conc"));
      b.setAttribute("aria-pressed", String(Math.abs(Math.log10(step) - Number(log)) < 0.02));
    });
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
        yaxis: axis({ title: "engagement (occupancy only for Kd/Ki rows)", range: [0, 1.02] }),
        xaxis: axis({ title: "" }),
      }
    );
    $("#legend-scorecard").innerHTML =
      `<span class="item" style="color:${c.insect}"><span class="swatch" style="background:${c.insect}"></span>insect target</span>` +
      `<span class="item" style="color:${c.vertebrate}"><span class="swatch" style="background:${c.vertebrate}"></span>vertebrate counterpart</span>` +
      `<span class="item">hatched / missing bar = no sourced value: not modelled, not zero</span>`;

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
          engagement(p.log10_ec50_ratio_vert_over_insect, 2),
          engagement(p.occupancy_difference, 3),
          `<span class="badge tier-${esc(p.evidence_tier)}">${esc(p.evidence_tier.replace(/_/g, " "))}</span>`,
        ];
      })
    );

    table(
      "tbl-occupancy",
      [
        { label: "receptor" },
        { label: "engagement", num: true },
        { label: "value (M)", num: true },
        { label: "type" },
        { label: "n", num: true },
        { label: "direction" },
        { label: "tier" },
      ],
      (occ.receptors || []).map((r) => {
        const m = receptorMeta(r.receptor);
        return [
          `<span class="swatch" style="background:${m.organism === "insect" ? c.insect : c.vertebrate}"></span>${esc(r.receptor)}`,
          engagement(r.engagement === undefined ? r.occupancy : r.engagement, 3),
          sci(r.param_value_M === undefined ? r.ec50_M : r.param_value_M),
          esc(r.param_type || ""),
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
    const receptors = Object.keys(pts[0].receptors || {}).filter(
      (k) => pts.some((p) => p.receptors[k] !== null && p.receptors[k] !== undefined)
    );
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
      yaxis: axis({ title: "engagement (occupancy only for Kd/Ki rows)", range: [0, 1.02] }),
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
    const fallbackT = isStatic() ? STATIC_T_MS : 500;
    const t_ms = Math.min(3000, Math.max(50, Number($("#f-tms").value) || fallbackT));
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

    // not-modelled receptors have no curve at all (null), so they are not drawn
    const receptors = Object.keys(out.occupancy || {}).filter((k) => out.occupancy[k]);
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
    persistDesign({ experiment: spec, experimentSettings: experimentSettings() });
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
    if (!out) {
      table(
        "tbl-experiment",
        [{ label: "compound" }, { label: "conc (M)", num: true }, { label: "n", num: true }],
        [],
        { empty: "Run the design to fill this grid." }
      );
      draw("plot-exp-heat", [], { height: 220 });
      $("#exp-status").textContent = "Preset loaded; run the design to calculate these settings.";
      return;
    }
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
  /* In the browser build every one of these shuffles is a circuit run on this
     machine, single-threaded, with the page frozen while it happens. Say so
     before starting instead of after. */
  function confirmLongRun(what, estimate) {
    if (!isStatic()) return true;
    const ok = window.confirm(
      `${what} runs entirely in this browser tab and takes about ${estimate}. ` +
        "The page will not respond until it finishes.\n\nRun it now?"
    );
    if (!ok) toast(`${what} cancelled`, "info");
    return ok;
  }

  async function runControls() {
    const d = design();
    if (!confirmLongRun("The null panel and the selectivity landscape", "1-3 minutes")) return null;
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
    if (!confirmLongRun(`A ${n}-shuffle ${mode} null distribution`, "10-60 seconds")) return;
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
  // 0. DASHBOARD - the interpretation layer
  // ==================================================================
  /* Two vocabularies, both from flylab/analysis/claims.py: CLASSIFICATIONS is
     what a reader sees on a chip, LABELS is what the machine-readable claim
     chain carries. Nothing on this page may print a verdict ("safe", "toxic",
     "87% confident"); it prints values, and says where each one came from. */
  const CHIPS = ["LITERATURE", "MODEL-DERIVED", "MODEL-ASSUMPTION", "PREDICTION", "NOT MODELLED"];
  const CHIP_NOTE = {
    LITERATURE: "a measured value from a cited paper (a potency, an affinity, or the MaleCNS reconstruction itself)",
    "MODEL-DERIVED": "computed from a literature value by a transformation FlyLab chose, such as the Hill engagement or a selectivity index",
    "MODEL-ASSUMPTION": "asserted by the model and never fitted to animal data: the gain rules, uniform expression, the engine constants",
    PREDICTION: "simulation output, including the predicted transmitter labels the signs rest on",
    "NOT MODELLED": "no sourced value exists; the row reports nothing rather than zero",
  };

  function chip(kind, provKey, label) {
    const k = CHIPS.indexOf(kind) >= 0 ? kind : "NOT MODELLED";
    const text = esc(label || k.toLowerCase());
    if (!provKey)
      return (
        `<button type="button" class="chip static" data-k="${esc(k)}" data-explain="${esc(k)}" ` +
        `title="${esc(CHIP_NOTE[k])} — click for what this label means">${text}</button>`
      );
    return (
      `<button type="button" class="chip" data-k="${esc(k)}" data-prov="${esc(provKey)}" ` +
      `title="${esc(CHIP_NOTE[k])} — click for provenance">${text}</button>`
    );
  }

  function provPut(key, record) {
    state.prov[key] = record;
    return key;
  }

  function provRow(label, value, cls) {
    if (value === null || value === undefined || value === "") return "";
    return `<dt>${esc(label)}</dt><dd class="${cls || ""}">${esc(value)}</dd>`;
  }

  function openProvenance(key) {
    const rec = state.prov[key];
    const drawer = $("#provenance");
    const body = $("#prov-body");
    const title = $("#prov-title");
    if (!drawer || !body) return;
    if (!rec) {
      body.innerHTML = '<p class="empty">No provenance recorded for this value.</p>';
    } else if (rec.kind === "receptor") {
      const r = rec.row || {};
      title.textContent = "Where this number came from";
      body.innerHTML =
        `<div class="prov-head"><span class="name">${esc(r.receptor)}</span>${chip(r.classification)}</div>` +
        "<dl>" +
        provRow("parameter type", r.param_type) +
        provRow("value", r.param_value_M === null || r.param_value_M === undefined ? "not modelled" : `${Number(r.param_value_M).toExponential(3)} M`) +
        provRow("Hill n", r.n) +
        provRow("direction", r.direction) +
        provRow("engagement at this dose", r.engagement === null || r.engagement === undefined ? "not modelled" : Number(r.engagement).toFixed(4)) +
        provRow("engagement model", r.engagement_model) +
        provRow("what that means", r.engagement_model_note) +
        provRow("evidence distance", r.evidence_distance_label || r.evidence_distance) +
        provRow("distance means", r.evidence_distance_note) +
        provRow("provenance warning", r.provenance_warning) +
        provRow("what the parameter is", r.param_type_note) +
        provRow("source relation", r.relation) +
        provRow("relation means", r.relation_note) +
        provRow("species / preparation", r.species) +
        provRow("evidence tier", r.evidence_tier) +
        provRow("why not modelled", r.not_modelled_reason) +
        provRow("source", r.source, "source") +
        (r.doi ? `<dt>DOI</dt><dd><a href="https://doi.org/${esc(r.doi)}" rel="noreferrer noopener" target="_blank">${esc(r.doi)}</a></dd>` : "") +
        (r.pmid ? `<dt>PMID</dt><dd><a href="https://pubmed.ncbi.nlm.nih.gov/${esc(r.pmid)}/" rel="noreferrer noopener" target="_blank">${esc(r.pmid)}</a></dd>` : "") +
        "</dl>";
    } else if (rec.kind === "claim") {
      const l = rec.link || {};
      title.textContent = "Claim chain link";
      body.innerHTML =
        `<div class="prov-head"><span class="name">${esc(l.step)}</span>${chip(l.classification)}` +
        `<span class="badge">${esc(l.label)}</span></div>` +
        `<p>${esc(l.statement)}</p>` +
        ((l.assumptions || []).length
          ? "<dt>assumptions it introduces</dt><dd><ul>" +
            l.assumptions.map((a) => `<li>${esc(a)}</li>`).join("") +
            "</ul></dd>"
          : "") +
        ((l.unknowns || []).length
          ? "<dt>what it leaves unknown</dt><dd><ul>" +
            l.unknowns.map((a) => `<li>${esc(a)}</li>`).join("") +
            "</ul></dd>"
          : "") +
        ((l.sources || []).length
          ? "<dt>sources</dt><dd class=\"source\">" +
            l.sources
              .map((sc) => `<p>${esc(sc.receptor ? sc.receptor + ": " : "")}${esc(sc.source || "")}</p>`)
              .join("") +
            "</dd>"
          : "") +
        `<dt>detail</dt><dd><pre class="json">${esc(JSON.stringify(l.detail || {}, null, 2))}</pre></dd>`;
    } else {
      title.textContent = rec.title || "Provenance";
      body.innerHTML =
        `<div class="prov-head"><span class="name">${esc(rec.title || "")}</span>${chip(rec.classification)}</div>` +
        `<p>${esc(rec.basis || rec.statement || "")}</p>` +
        `<dt>detail</dt><dd><pre class="json">${esc(JSON.stringify(rec.detail || {}, null, 2))}</pre></dd>`;
    }
    drawer.dataset.open = "true";
    drawer.setAttribute("aria-hidden", "false");
  }

  function closeProvenance() {
    const drawer = $("#provenance");
    if (!drawer) return;
    drawer.dataset.open = "false";
    drawer.setAttribute("aria-hidden", "true");
  }

  // ---- data --------------------------------------------------------
  /* The browser build asks for FAST_N shuffles and says so in the caption:
     a dependence verdict at n=20 resolves p no finer than 1/(n+1). */
  const BROWSER_DEPENDENCE_N = 20;

  async function runDashboard() {
    const d = design();
    const out = await api(
      `/api/dashboard?compound=${encodeURIComponent(d.compound)}&conc_M=${d.conc_M}` +
        `&graph=${encodeURIComponent(d.graph)}&n=${BROWSER_DEPENDENCE_N}` +
        (d.genotype ? `&genotype=${encodeURIComponent(d.genotype)}` : "")
    );
    state.cache.dashboard = out;
    if (out.circuit && out.circuit.notebook) setNotebook(out.circuit.notebook, "dashboard · subgraph");
    renderDashboard();
    return (out.runtime_s || 0) * 1000;
  }

  function indexProvenance(out) {
    state.prov = state.prov || {};
    (out.evidence || []).forEach((r) => provPut("receptor:" + r.receptor, { kind: "receptor", row: r }));
    ((out.claims || {}).chain || []).forEach((l) => provPut("claim:" + l.step, { kind: "claim", link: l }));
    ((out.trust || {}).rows || []).forEach((r) =>
      provPut("trust:" + r.area, {
        kind: "trust",
        title: r.area,
        classification: r.classification,
        basis: r.basis,
        detail: r.detail,
      })
    );
  }

  // ---- render ------------------------------------------------------
  function renderDashboard() {
    const out = state.cache.dashboard;
    if (!out) return;
    indexProvenance(out);
    renderOverview(out);
    renderTiles(out);
    renderWhy(out, "dash-why");
    renderLadder(out, "plot-dash-ladder", "legend-dash-ladder", 300);
    renderLadderTable(out);
    renderSelectivityPanel(out);
    renderCircuitConsequence(out);
    renderDependence(out);
    renderTrust(out);
    if (state.cache.compare) renderCompare();
  }

  function renderOverview(out) {
    const host = $("#dash-overview");
    if (!host) return;
    const c = out.compound || {};
    const cov = out.coverage || {};
    const conc = Number(out.concentration_M);
    const dose =
      conc >= 1e-6 ? `${(conc * 1e6).toPrecision(3)} µM` : conc >= 1e-9 ? `${(conc * 1e9).toPrecision(3)} nM` : `${conc.toExponential(2)} M`;
    const target = c.target_receptor
      ? `${esc(String(c.target_receptor).replace("insect_", ""))} <span class="muted">(${esc(c.mode || "—")})</span>`
      : '<span class="muted">no sourced insect target</span>';
    const warn = (out.warnings || []).slice(0, 3);
    host.className = "overview";
    host.innerHTML =
      `<p class="hero">${esc(c.name || c.key)} <span class="dose">· ${esc(dose)}</span></p>` +
      `<div class="facts">` +
      `<span>class <b>${esc(c.class || "unclassified")}</b></span>` +
      `<span>target <b>${target}</b></span>` +
      (c.cas ? `<span>CAS <b>${esc(c.cas)}</b></span>` : "") +
      `<span>graph <b>${esc(out.graph)}</b></span>` +
      `</div>` +
      `<div class="coverage">Evidence coverage: <b>${cov.n_sourced}</b> of <b>${cov.n_rows}</b> receptor rows are literature-sourced, ` +
      `<b>${cov.n_not_modelled}</b> are not modelled (no value exists — they are not zero). ` +
      `${chip("LITERATURE", "claim:parameter_source", "sourced")} ${chip("NOT MODELLED", null, "not modelled")}</div>` +
      (warn.length ? warn.map((w) => `<div class="notice">${esc(w)}</div>`).join("") : "");
  }

  function tileHtml(label, value, unit, kind, provKey, note) {
    return (
      `<div class="tile"><div class="k">${esc(label)}</div>` +
      `<div class="v big">${value}</div>` +
      `<div class="u">${esc(unit || "")}</div>` +
      (note ? `<div class="note">${esc(note)}</div>` : "") +
      `<div class="chiprow">${chip(kind, provKey)}</div></div>`
    );
  }

  function renderTiles(out) {
    const host = $("#dash-tiles");
    if (!host) return;
    const h = out.headline || {};
    const ie = h.insect_engagement || {};
    const ve = h.vertebrate_engagement || {};
    const rs = h.receptor_selectivity || {};
    const et = h.evidence_tier || {};
    host.innerHTML =
      tileHtml(
        "Insect engagement",
        engagement(ie.value, 3),
        ie.receptor ? String(ie.receptor).replace("insect_", "") + " · " + (ie.param_type || "") : "no sourced target",
        ie.classification,
        ie.receptor ? "receptor:" + ie.receptor : null
      ) +
      tileHtml(
        "Vertebrate engagement",
        engagement(ve.value, 3),
        ve.receptor ? String(ve.receptor).replace("vertebrate_", "") + " · " + (ve.param_type || "") : "no sourced target",
        ve.classification,
        ve.receptor ? "receptor:" + ve.receptor : null
      ) +
      tileHtml(
        "Receptor selectivity",
        rs.ratio_vert_over_insect === null || rs.ratio_vert_over_insect === undefined
          ? NOT_MODELLED
          : Number(rs.ratio_vert_over_insect).toPrecision(3) + "&times;",
        rs.pair ? rs.pair + " pair · potency ratio" : "no comparable pair",
        rs.classification,
        "claim:engagement_transformation",
        "ratio of potencies — not a safety margin"
      ) +
      tileHtml(
        "Evidence tier",
        esc(String(et.value || "—").replace(/_/g, " ")),
        `${et.n_sourced} sourced / ${et.n_not_modelled} not modelled`,
        et.classification,
        "claim:parameter_source"
      );
  }

  function renderWhy(out, hostId) {
    const host = document.getElementById(hostId);
    if (!host) return;
    const why = out.why || {};
    const inputs = why.inputs || {};
    host.innerHTML =
      `<h3>Why this happened</h3>` +
      `<p>${esc(why.text || "—")}</p>` +
      `<div class="src">Generated from this run — saturation ${esc(inputs.saturation || "—")}, ` +
      `dependence class ${esc(inputs.dependence_class || "not run")}, ` +
      `dominant term ${esc(inputs.dominant_uncertainty || "—")}. ${chip("MODEL-DERIVED", "claim:mechanism_rule")}</div>`;
  }

  // ---- block 4: one concentration axis ------------------------------
  function ladderTraces(out) {
    const lad = out.ladder;
    if (!lad) return null;
    const c = colors();
    const x = (lad.concs_M || []).map((v) => Math.log10(v));
    const series = [
      { name: "insect engagement", y: lad.insect_engagement, color: c.insect, dash: "solid" },
      { name: "vertebrate engagement", y: lad.vertebrate_engagement, color: c.vertebrate, dash: "dash" },
      { name: "circuit response", y: lad.circuit_response, color: cssVar("--nt-glutamate"), dash: "solid" },
    ];
    const traces = series
      .filter((s) => (s.y || []).some((v) => v !== null && v !== undefined))
      .map((s) => ({
        type: "scatter",
        mode: "lines",
        name: s.name,
        x: x,
        y: s.y,
        line: { color: s.color, width: 2, shape: "spline", smoothing: 0.4, dash: s.dash },
        hovertemplate: `${esc(s.name)}<br>10^%{x:.1f} M · %{y:.3f}<extra></extra>`,
      }));
    return { traces: traces, x: x, series: series, lad: lad };
  }

  function renderLadder(out, plotId, legendId, height) {
    const built = ladderTraces(out);
    if (!built) return;
    const { traces, x, series, lad } = built;
    const c = colors();
    const shapes = [];
    const annotations = [];
    /* Three vertical rules can land within a decade of each other, so their
       labels are stepped down the plot instead of stacked on one line. Series
       identity comes from the legend below the chart: the three curves converge
       at the right edge, where end-labels would detach from their lines. */
    const rule = (value, label, color, slot) => {
      if (!value || !Number.isFinite(Math.log10(value))) return;
      const lx = Math.log10(value);
      if (lx < x[0] || lx > x[x.length - 1]) return;
      shapes.push({ type: "line", x0: lx, x1: lx, yref: "paper", y0: 0, y1: 1, line: { color: color, width: 2 } });
      const right = lx > (x[0] + x[x.length - 1]) / 2;
      annotations.push({
        x: lx,
        y: [1.0, 0.9, 0.8][slot] || 0.7,
        text: (right ? label + " " : " " + label),
        showarrow: false,
        yanchor: "top",
        xanchor: right ? "right" : "left",
        font: { color: color, size: 10.5 },
      });
    };
    rule(lad.circuit_threshold_M, "circuit 50%", cssVar("--nt-glutamate"), 0);
    rule(lad.vertebrate_threshold_M, "vertebrate 20%", c.vertebrate, 1);
    rule(lad.current_conc_M, "current dose", c.muted, 2);
    draw(plotId, traces, {
      showlegend: true,
      height: height || 320,
      shapes: shapes,
      annotations: annotations,
      margin: { t: 34, r: 18, b: 44, l: 56 },
      xaxis: axis({ title: "log10 concentration (M)" }),
      yaxis: axis({ title: "fraction (0–1)", range: [0, 1.04] }),
    });
    const host = document.getElementById(legendId);
    if (host) {
      host.innerHTML =
        `<span class="item" style="color:${c.insect}"><span class="solid"></span>insect engagement</span>` +
        `<span class="item" style="color:${c.vertebrate}"><span class="dash"></span>vertebrate engagement</span>` +
        `<span class="item" style="color:${cssVar("--nt-glutamate")}"><span class="solid"></span>circuit response = |treated − vehicle| / vehicle</span>` +
        `<span class="item">readout ${esc(lad.readout || "")} on ${esc(out.graph)} · one axis, no second scale</span>`;
    }
  }

  function renderLadderTable(out) {
    const lad = out.ladder;
    if (!lad) return;
    table(
      "tbl-dash-ladder",
      [
        { label: "conc (M)", num: true },
        { label: "insect engagement", num: true },
        { label: "vertebrate engagement", num: true },
        { label: "circuit response", num: true },
        { label: "treated (Hz)", num: true },
      ],
      (lad.concs_M || []).map((conc, i) => [
        sci(conc, 1),
        engagement(lad.insect_engagement[i], 3),
        engagement(lad.vertebrate_engagement[i], 3),
        engagement(lad.circuit_response[i], 3),
        num(lad.circuit_treated_hz[i], 3),
      ]),
      { empty: "No ladder for this compound." }
    );
  }

  // ---- block 2: selectivity ----------------------------------------
  function meter(name, value, color, max) {
    const v = Number(value);
    const pct = Number.isFinite(v) ? Math.max(0, Math.min(1, v / (max || 1))) * 100 : 0;
    return (
      `<div class="meter"><span class="name">${esc(name)}</span>` +
      `<span class="track"><span class="fill" style="left:0;width:${pct.toFixed(1)}%;background:${color}"></span></span>` +
      `<span class="val">${Number.isFinite(v) ? v.toFixed(3) : "—"}</span></div>`
    );
  }

  function renderSelectivityPanel(out) {
    const sel = out.selectivity || {};
    const c = colors();
    const host = $("#meters-dash-sel");
    if (host) {
      host.innerHTML =
        meter(
          String(sel.insect_receptor || "insect").replace("insect_", ""),
          sel.insect_engagement,
          c.insect
        ) +
        meter(
          String(sel.vertebrate_receptor || "vertebrate").replace("vertebrate_", ""),
          sel.vertebrate_engagement,
          c.vertebrate
        );
    }
    const facts = $("#dash-sel-facts");
    const lim = sel.limiting_vertebrate_receptor || null;
    if (facts) {
      facts.innerHTML =
        `<div class="facts-list">` +
        `<p><b>${sel.ratio_vert_over_insect === null || sel.ratio_vert_over_insect === undefined ? "—" : Number(sel.ratio_vert_over_insect).toPrecision(3) + "×"}</b> ` +
        `insect/vertebrate potency ratio on the ${esc(sel.pair || "—")} pair ${chip("MODEL-DERIVED", "claim:engagement_transformation")}</p>` +
        `<p><b>${engagement(sel.engagement_difference, 3)}</b> engagement difference at this dose ` +
        `(insect ${engagement(sel.insect_engagement, 3)} − vertebrate ${engagement(sel.vertebrate_engagement, 3)})</p>` +
        `<p><b>${sci(sel.vertebrate_limit_conc_M, 2)} M</b> is where vertebrate engagement first reaches ` +
        `${Math.round((sel.occ_limit || 0.2) * 100)}%` +
        (lim
          ? `, and <b>${esc(String(lim.receptor).replace("vertebrate_", ""))}</b> is the receptor that becomes limiting first ` +
            chip(lim.classification, "receptor:" + lim.receptor)
          : "") +
        `</p>` +
        (sel.dose_over_vertebrate_limit
          ? `<p>The current dose is <b>${Number(sel.dose_over_vertebrate_limit).toPrecision(3)}×</b> that concentration.</p>`
          : "") +
        `<p class="hint">${esc(sel.statement || "")}</p>` +
        `</div>`;
    }
    table(
      "tbl-dash-vert",
      [
        { label: "vertebrate receptor" },
        { label: "engagement here", num: true },
        { label: "sourced value", num: true },
        { label: "20% at (M)", num: true },
        { label: "evidence" },
      ],
      (sel.vertebrate_rows || []).map((r) => [
        `<span class="swatch" style="background:${c.vertebrate}"></span>${esc(String(r.receptor).replace("vertebrate_", ""))}`,
        engagement(r.engagement_at_dose, 3),
        r.param_value_M === null || r.param_value_M === undefined
          ? NOT_MODELLED
          : `${esc(r.param_type)} ${sci(r.param_value_M, 2)}`,
        sci(r.conc_at_limit_M, 2),
        chip(r.classification, "receptor:" + r.receptor),
      ]),
      { empty: "This compound lists no vertebrate receptor." }
    );
  }

  // ---- block 3: circuit consequence ---------------------------------
  function renderCircuitConsequence(out) {
    const cir = out.circuit || {};
    const rows = cir.readouts || {};
    const keys = Object.keys(rows);
    const c = colors();
    const labels = keys.map((k) => rows[k].label);
    // dumbbell: before -> after per item, one hue in two shades
    const traces = [
      {
        type: "scatter",
        mode: "lines",
        x: keys.flatMap((k) => [rows[k].vehicle, rows[k].treated, null]),
        y: keys.flatMap((k) => [rows[k].label, rows[k].label, null]),
        line: { color: c.grid, width: 2 },
        hoverinfo: "skip",
        showlegend: false,
      },
      {
        type: "scatter",
        mode: "markers",
        name: "vehicle",
        x: keys.map((k) => rows[k].vehicle),
        y: labels,
        marker: { color: cssVar("--seq-250"), size: 12, line: { width: 2, color: c.surface } },
        hovertemplate: "vehicle %{x:.3f} Hz<extra></extra>",
      },
      {
        type: "scatter",
        mode: "markers+text",
        name: "treated",
        x: keys.map((k) => rows[k].treated),
        y: labels,
        text: keys.map((k) => (rows[k].percent === null ? "" : `${rows[k].percent > 0 ? "+" : ""}${rows[k].percent.toFixed(1)}%`)),
        textposition: "middle right",
        textfont: { color: c.text2, size: 10.5 },
        marker: { color: cssVar("--seq-550"), size: 12, line: { width: 2, color: c.surface } },
        hovertemplate: "treated %{x:.3f} Hz<extra></extra>",
      },
    ];
    draw("plot-dash-circuit", traces, {
      showlegend: true,
      height: 250,
      margin: { t: 34, r: 60, b: 44, l: 120 },
      xaxis: axis({ title: "firing rate (Hz, simulated)" }),
      yaxis: axis({ title: "", automargin: true }),
    });
    const legend = $("#legend-dash-circuit");
    if (legend) {
      legend.innerHTML =
        `<span class="item" style="color:${cssVar("--seq-250")}"><span class="swatch" style="background:${cssVar("--seq-250")}"></span>vehicle</span>` +
        `<span class="item" style="color:${cssVar("--seq-550")}"><span class="swatch" style="background:${cssVar("--seq-550")}"></span>treated</span>` +
        `<span class="item">${chip("PREDICTION", "claim:readout")} every rate here is model output</span>`;
    }
    table(
      "tbl-dash-circuit",
      [
        { label: "readout" },
        { label: "vehicle (Hz)", num: true },
        { label: "treated (Hz)", num: true },
        { label: "Δ (Hz)", num: true },
        { label: "change", num: true },
        { label: "direction" },
      ],
      keys.map((k) => [
        esc(rows[k].label),
        num(rows[k].vehicle, 3),
        num(rows[k].treated, 3),
        num(rows[k].delta, 3),
        rows[k].percent === null ? "—" : `${rows[k].percent > 0 ? "+" : ""}${num(rows[k].percent, 1)}%`,
        esc(rows[k].direction),
      ]),
      { empty: "No circuit readouts." }
    );
  }

  function renderDependence(out) {
    const dep = out.dependence;
    const host = $("#dash-dep-verdict");
    if (!dep) {
      if (host) host.innerHTML = '<div class="empty">Dependence not requested for this run.</div>';
      table("tbl-dash-dep", [{ label: "" }], [], { empty: "Not run." });
      return;
    }
    if (host) {
      host.innerHTML =
        `<div class="verdict" data-class="${esc(dep.class)}"><span class="mark"></span>` +
        `<span>${esc(dep.verdict)}</span>${chip("MODEL-DERIVED", "claim:malecns_edges")}` +
        `<span class="why">${esc(dep.reason || "")}</span></div>` +
        `<div class="hint">${dep.n} shuffles per mode, so the empirical p resolves no finer than ` +
        `${num(dep.p_resolution, 4)}. Browser default is the fast setting; raise n for a quotable p. ` +
        `${esc(dep.note || "")}</div>`;
    }
    table(
      "tbl-dash-dep",
      [
        { label: "null model" },
        { label: "information it keeps" },
        { label: "p (two-sided)", num: true },
        { label: "resolution", num: true },
        { label: "verdict" },
      ],
      (dep.modes || []).map((m) => [
        `<span class="mono">${esc(m.mode)}</span>`,
        esc(m.information_kept || ""),
        (m.at_resolution_floor ? "≤ " : "") + num(m.p_two_sided, 4),
        num(m.p_resolution, 4),
        depVerdictBadge(m),
      ]),
      { empty: "No null modes ran." }
    );
  }

  /* The three verdicts of flylab/analysis/dependence.py mode_verdict(), never
     two. A failure to reject is "not distinguishable"; only a gap from the
     null median below the prespecified margin is equivalence. The older label
     here said "shuffle reproduces it" for every non-rejection, which is the
     one sentence NOVELTY.md forbids. If the payload carries the per-mode
     verdict we print it; if it carries only beats_null we print the weaker of
     the two answers it can support, which is the correct one. */
  const DEP_VERDICT_LABEL = {
    distinguishable: ["real effect stands out", "tier-literature_order",
      "the permutation test rejects: the real graph and this graph model give different drug effects"],
    equivalent_within_tolerance: ["same within margin", "tier-measured_fit",
      "the test did not reject AND the gap from the null median is below the prespecified margin"],
    indeterminate: ["not distinguishable", "",
      "a failure to reject, and the gap from the null median is not below the prespecified margin: consistent with a difference this run cannot resolve"],
  };
  function depVerdictBadge(m) {
    const v =
      m.verdict ||
      (m.beats_null ? "distinguishable" : m.equivalent_within_tolerance ? "equivalent_within_tolerance" : null);
    if (v && DEP_VERDICT_LABEL[v]) {
      const [text, cls, why] = DEP_VERDICT_LABEL[v];
      return `<span class="badge ${cls}" title="${esc(why)}">${esc(text)}</span>`;
    }
    // beats_null alone: "not distinguishable" is everything it licenses
    return (
      '<span class="badge" title="a failure to reject. Whether it also reaches equivalence within the ' +
      'prespecified margin is stated in the verdict line above; a non-rejection on its own is not ' +
      'evidence that the shuffle gives the same effect.">not distinguishable</span>'
    );
  }

  // ---- block 8: what can I trust? ------------------------------------
  function renderTrust(out) {
    const trust = out.trust || {};
    table(
      "tbl-dash-trust",
      [{ label: "layer" }, { label: "status" }, { label: "on what basis" }, { label: "classification" }],
      (trust.rows || []).map((r) => [
        `<b>${esc(r.area)}</b>`,
        esc(r.status),
        esc(r.basis),
        chip(r.classification, "trust:" + r.area),
      ]),
      { empty: "Not run yet." }
    );
    const note = $("#trust-note");
    if (note) note.textContent = trust.note || "";
  }

  // ---- block 9: compare compounds ------------------------------------
  async function estimateCompare() {
    const picked = $$("#f-cmp-compounds option")
      .filter((o) => o.selected)
      .map((o) => o.value);
    const host = $("#cmp-estimate");
    if (!host) return picked;
    if (!picked.length) {
      host.textContent = "pick two or more compounds";
      return picked;
    }
    try {
      const est = await post("/api/compare", {
        compounds: picked,
        conc_M: design().conc_M,
        include_dependence: $("#f-cmp-dependence").checked,
        n: BROWSER_DEPENDENCE_N,
        estimate_only: true,
      });
      const s = (est.runtime_estimate || {}).estimate_s;
      host.textContent = `estimated ${s} s for ${picked.length} compounds`;
    } catch (err) {
      host.textContent = "";
    }
    return picked;
  }

  async function runCompare() {
    const picked = await estimateCompare();
    if (picked.length < 1) throw new Error("pick at least one compound");
    const out = await post("/api/compare", {
      compounds: picked,
      conc_M: design().conc_M,
      include_dependence: $("#f-cmp-dependence").checked,
      n: BROWSER_DEPENDENCE_N,
    });
    state.cache.compare = out;
    renderCompare();
    return (out.runtime_s || 0) * 1000;
  }

  function renderCompare() {
    const out = state.cache.compare;
    if (!out) return;
    const c = colors();
    const rows = out.rows || [];
    const scored = rows.filter(
      (r) => r.receptor_si_log10 !== null && r.circuit_si_log10 !== null
    );
    const lo = Math.min(0, ...scored.map((r) => Math.min(r.receptor_si_log10, r.circuit_si_log10)));
    const hi = Math.max(1, ...scored.map((r) => Math.max(r.receptor_si_log10, r.circuit_si_log10)));
    draw(
      "plot-compare",
      [
        {
          type: "scatter",
          mode: "markers+text",
          x: scored.map((r) => r.receptor_si_log10),
          y: scored.map((r) => r.circuit_si_log10),
          text: scored.map((r) => r.name),
          textposition: "top center",
          textfont: { color: c.text2, size: 10.5 },
          cliponaxis: false,
          marker: { color: cssVar("--seq-550"), size: 12, line: { width: 2, color: c.surface } },
          hovertemplate: "%{text}<br>receptor SI %{x:.2f}<br>circuit SI %{y:.2f}<extra></extra>",
        },
      ],
      {
        height: 320,
        margin: { t: 40, r: 30, b: 44, l: 56 },
        xaxis: axis({ title: "receptor selectivity (log10)", range: [lo - 0.5, hi + 0.5] }),
        yaxis: axis({ title: "circuit selectivity (log10)", range: [lo - 0.5, hi + 0.9] }),
        shapes: [
          {
            type: "line",
            x0: lo - 0.5,
            y0: lo - 0.5,
            x1: hi + 0.5,
            y1: hi + 0.5,
            line: { color: c.base, width: 1 },
          },
        ],
        annotations: [
          {
            x: hi + 0.5,
            y: hi + 0.5,
            text: "equal selectivity",
            showarrow: false,
            xanchor: "right",
            yanchor: "bottom",
            font: { color: c.muted, size: 10.5 },
          },
        ],
      }
    );
    const findings = $("#cmp-findings");
    if (findings) {
      const missing = rows.filter((r) => r.circuit_si_log10 === null).map((r) => r.name);
      findings.innerHTML =
        (out.findings || []).map((f) => `<div class="notice info">${esc(f)}</div>`).join("") +
        (missing.length
          ? `<div class="notice">${esc(missing.join(", "))} ${missing.length > 1 ? "are" : "is"} ` +
            "not on the plot: the circuit never reaches a 50% change on this cut, so no circuit " +
            "selectivity index exists. That is a result, not a missing value.</div>"
          : "");
    }
    table(
      "tbl-compare",
      [
        { label: "compound" },
        { label: "target" },
        { label: "insect eng.", num: true },
        { label: "vert eng.", num: true },
        { label: "receptor SI", num: true },
        { label: "circuit SI", num: true },
        { label: "SI gap", num: true },
        { label: "circuit Δ", num: true },
        { label: "topology dependence" },
        { label: "evidence" },
      ],
      rows.map((r) => [
        `<b>${esc(r.name)}</b>`,
        esc(String(r.target_receptor || "—").replace("insect_", "")),
        engagement(r.insect_engagement, 3),
        engagement(r.vertebrate_engagement, 3),
        r.receptor_si_log10 === null ? NOT_MODELLED : num(r.receptor_si_log10, 2),
        r.circuit_si_log10 === null
          ? '<span class="muted">50% never reached</span>'
          : num(r.circuit_si_log10, 2),
        r.si_gap_circuit_minus_receptor === null ? "—" : num(r.si_gap_circuit_minus_receptor, 2),
        r.circuit_delta_percent === null
          ? "—"
          : `${r.circuit_delta_percent > 0 ? "+" : ""}${num(r.circuit_delta_percent, 1)}%`,
        esc(r.topology_dependence || "not run"),
        chip(
          r.evidence_tier === "literature_order" ? "LITERATURE" : "MODEL-ASSUMPTION",
          null,
          `${r.n_sourced}/${r.n_rows} sourced`
        ),
      ]),
      { empty: "Pick compounds and press Run compare." }
    );
  }

  // ==================================================================
  // 0b. DOSE-RESPONSE
  // ==================================================================
  async function runDose() {
    /* The dose panel is a second view of the dashboard payload: one call feeds
       both, so opening either after the other costs nothing. */
    if (!state.cache.dashboard) {
      await runDashboard();
      state.loaded.dashboard = true;
    }
    renderDose();
    return null;
  }

  function renderDose() {
    const out = state.cache.dashboard;
    if (!out) return;
    renderLadder(out, "plot-dose-ladder", "legend-dose-ladder", 420);
    const inputs = (out.why || {}).inputs || {};
    const span = inputs.two_fold_engagement_span;
    const tiles = $("#tiles-saturation");
    if (tiles) {
      tiles.innerHTML =
        tile("insect engagement", engagement(inputs.insect_engagement, 3), inputs.saturation || "") +
        tile(
          "two-fold parameter error",
          span === null || span === undefined ? "—" : "±" + num(span, 3),
          "moves engagement by this much"
        ) +
        tile("dominant term", esc(inputs.potency_is_inert ? "gain rule" : "potency + gain rule"), "at this dose") +
        tile("dependence", esc(inputs.dependence_class || "not run"), "permutation nulls");
    }
    const text = $("#saturation-text");
    if (text) {
      const receptor = String(inputs.insect_receptor || "").replace("insect_", "");
      text.innerHTML =
        `<p>${
          inputs.potency_is_inert
            ? `Halving or doubling the cited ${esc(inputs.param_type || "potency")} for ${esc(receptor)} moves engagement by ${num(span, 3)} — ` +
              "effectively nothing. A saturated receptor cannot report a potency change, so at this dose the model's answer is set by " +
              "the occupancy-to-gain rule, not by the library value. Move the slider left until the curve leaves its plateau to see the potency matter again."
            : `A two-fold error in the cited ${esc(inputs.param_type || "potency")} for ${esc(receptor)} moves engagement by ${num(span, 3)}, ` +
              "so this dose sits on the informative part of the curve and the library value still matters."
        }</p>` +
        `<div class="src">${chip("MODEL-DERIVED", "claim:engagement_transformation")} derived from this run's Hill parameters.</div>`;
    }
    const lad = out.ladder || {};
    const c = colors();
    const x = (lad.concs_M || []).map((v) => Math.log10(v));
    /* Colour follows the receptor FAMILY, so nAChR is the same hue everywhere in
       the bench. Two receptors of one family in one organism would then share a
       colour, so the line style is the second channel: insect solid then dotted,
       vertebrate dashed then dash-dot. The legend and the table carry the names. */
    const seen = {};
    const traces = (lad.per_receptor || [])
      .filter((r) => (r.values || []).some((v) => v !== null && v !== undefined))
      .map((r) => {
        const m = receptorMeta(r.receptor);
        const slot = m.family + ":" + r.organism;
        const i = seen[slot] === undefined ? (seen[slot] = 0) : (seen[slot] += 1);
        const dashes = r.organism === "insect" ? ["solid", "dot", "longdash"] : ["dash", "dashdot", "longdashdot"];
        return {
          type: "scatter",
          mode: "lines",
          name: r.receptor,
          x: x,
          y: r.values,
          line: {
            color: familyColor(m.family),
            width: 2,
            dash: dashes[Math.min(i, dashes.length - 1)],
          },
          hovertemplate: `${esc(r.receptor)}<br>10^%{x:.1f} M · %{y:.3f}<extra></extra>`,
        };
      });
    draw("plot-dose-receptors", traces, {
      showlegend: true,
      height: 340,
      margin: { t: 34, r: 18, b: 44, l: 56 },
      xaxis: axis({ title: "log10 concentration (M)" }),
      yaxis: axis({ title: "engagement (0–1)", range: [0, 1.04] }),
    });
    table(
      "tbl-dose-receptors",
      [
        { label: "receptor" },
        { label: "organism" },
        { label: "param" },
        { label: "value (M)", num: true },
        { label: "evidence" },
      ],
      (out.evidence || []).map((r) => [
        `<span class="swatch" style="background:${r.organism === "insect" ? c.insect : c.vertebrate}"></span>${esc(r.receptor)}`,
        esc(r.organism),
        esc(r.param_type || ""),
        sci(r.param_value_M, 2),
        chip(r.classification, "receptor:" + r.receptor),
      ]),
      { empty: "No rows." }
    );
  }

  // ==================================================================
  // 0c. GENOTYPE
  // ==================================================================
  async function runGenotype() {
    const d = design();
    const out = await post("/api/genotype/panel", { compound: d.compound, conc_M: d.conc_M });
    state.cache.genotype = { panel: out, other: state.cache.genotype && state.cache.genotype.other };
    renderGenotype();
    return null;
  }

  async function runGenotypeOther() {
    const other = $("#f-geno-other").value;
    if (!other) return;
    const out = await post("/api/genotype/panel", { compound: other, conc_M: design().conc_M });
    state.cache.genotype = Object.assign({}, state.cache.genotype, { other: out });
    renderGenotype();
  }

  function genotypeRows(panel) {
    return (panel && panel.rows) || [];
  }

  function renderGenotype() {
    const data = state.cache.genotype;
    if (!data || !data.panel) return;
    const rows = genotypeRows(data.panel);
    const c = colors();
    const labels = rows.map((r) => (r.allele === "wild type" ? "wild type" : `${r.gene || ""} ${r.allele}`.trim()));
    notices("genotype-notices", data.panel.warnings || [], "");
    draw(
      "plot-geno-occ",
      [
        {
          type: "bar",
          x: labels,
          y: rows.map((r) => r.occupancy),
          text: rows.map((r) => (r.occupancy === null ? "" : Number(r.occupancy).toFixed(3))),
          textposition: "outside",
          textfont: { color: c.text2, size: 10.5 },
          width: rows.length > 3 ? 0.45 : 0.2,
          marker: { color: cssVar("--seq-550"), cornerradius: 4, line: { width: 2, color: c.surface } },
          hovertemplate: "%{x}<br>engagement %{y:.3f}<extra></extra>",
        },
      ],
      {
        height: 280,
        bargap: 0.55,
        margin: { t: 26, r: 18, b: 60, l: 56 },
        xaxis: axis({ title: "" }),
        yaxis: axis({ title: "engagement at this dose (0–1)", range: [0, 1.12] }),
      }
    );
    const mn9 = rows.map((r) => (r.readouts || {}).mn9_hz);
    draw(
      "plot-geno-circuit",
      [
        {
          type: "bar",
          x: labels,
          y: mn9,
          text: mn9.map((v) => (v === null || v === undefined ? "" : Number(v).toFixed(1))),
          textposition: "outside",
          textfont: { color: c.text2, size: 10.5 },
          width: rows.length > 3 ? 0.45 : 0.2,
          marker: { color: cssVar("--seq-400"), cornerradius: 4, line: { width: 2, color: c.surface } },
          hovertemplate: "%{x}<br>MN9 %{y:.2f} Hz<extra></extra>",
        },
      ],
      {
        height: 280,
        bargap: 0.55,
        margin: { t: 26, r: 18, b: 60, l: 56 },
        xaxis: axis({ title: "" }),
        yaxis: axis({
          title: "MN9 rate (Hz, simulated)",
          range: [0, Math.max(1, ...mn9.filter((v) => Number.isFinite(v))) * 1.18],
        }),
      }
    );
    const head = [
      { label: "genotype" },
      { label: "gene" },
      { label: "fold shift", num: true },
      { label: "value (M)", num: true },
      { label: "engagement", num: true },
      { label: "MN9 (Hz)", num: true },
      { label: "mean (Hz)", num: true },
    ];
    const body = (panel) =>
      genotypeRows(panel).map((r) => [
        esc(r.allele === "wild type" ? "wild type" : r.allele),
        esc(r.gene || "—"),
        r.fold_shift === null || r.fold_shift === undefined ? "—" : `${num(r.fold_shift, 1)}×`,
        sci(r.ec50_M, 2),
        engagement(r.occupancy, 3),
        num((r.readouts || {}).mn9_hz, 2),
        num((r.readouts || {}).mean_hz, 3),
      ]);
    table("tbl-genotype", head, body(data.panel), { empty: "No allele in the library touches this compound." });
    table("tbl-genotype-other", head, body(data.other), {
      empty: "Pick a second compound and press Run to see the same alleles applied to it.",
    });
  }

  // ==================================================================
  // 0d. MIXTURES
  // ==================================================================
  async function runMixtures() {
    const a = $("#f-mix-a").value;
    const b = $("#f-mix-b").value;
    const ca = Number($("#f-mix-a-conc").value);
    const cb = Number($("#f-mix-b-conc").value);
    const model = $("#f-mix-model").value;
    const mix = await post("/api/mixture", {
      components: [
        { compound: a, conc_M: Number.isFinite(ca) ? ca : 0 },
        { compound: b, conc_M: Number.isFinite(cb) ? cb : 0 },
      ],
      assay: "subgraph",
      model: model,
    });
    let iso = null;
    try {
      iso = await post("/api/mixture/isobologram", { compound_a: a, compound_b: b });
    } catch (err) {
      iso = null;
    }
    state.cache.mixtures = { mix: mix, iso: iso, a: a, b: b };
    setNotebook(mix, `mixture · ${a} + ${b}`);
    renderMixtures();
    return null;
  }

  function renderMixtures() {
    const data = state.cache.mixtures;
    if (!data) return;
    const { mix, iso } = data;
    const c = colors();
    const r = mix.readouts || {};
    const syn = mix.synergy || {};
    const singles = r.singles || [];
    const observed = r.effect_fraction;
    const expected = expectedFromWhy(syn.why);
    notices("mixture-notices", mix.warnings || [], "");
    const tiles = $("#tiles-mixture");
    if (tiles) {
      tiles.innerHTML =
        tile("observed effect", engagement(observed, 3), "|treated − vehicle| / vehicle") +
        tile(`${esc(syn.model || "bliss")} expectation`, expected === null ? "—" : num(expected, 3), "null model") +
        tile("interaction", esc(syn.verdict || "—"), "classification") +
        singles
          .map((s) =>
            tile(
              esc(s.compound),
              engagement(s.effect_fraction === undefined ? s.relative_change : s.effect_fraction, 3),
              sci(s.conc_M, 1) + " M alone"
            )
          )
          .join("");
    }
    const labels = singles.map((s) => s.compound).concat(["combination"]);
    const values = singles
      .map((s) => (s.effect_fraction === undefined ? s.relative_change : s.effect_fraction))
      .concat([observed]);
    const shapes = [];
    const annotations = [];
    if (expected !== null) {
      shapes.push({
        type: "line",
        xref: "paper",
        x0: 0,
        x1: 1,
        y0: expected,
        y1: expected,
        line: { color: c.base, width: 1 },
      });

    }
    draw(
      "plot-mixture",
      [
        {
          type: "bar",
          x: labels,
          y: values,
          text: values.map((v) => (v === null || v === undefined ? "" : Number(v).toFixed(3))),
          textposition: "outside",
          textfont: { color: c.text2, size: 10.5 },
          width: labels.length > 3 ? 0.45 : 0.2,
          marker: { color: cssVar("--seq-550"), cornerradius: 4, line: { width: 2, color: c.surface } },
          hovertemplate: "%{x}<br>effect %{y:.3f}<extra></extra>",
        },
      ],
      {
        height: 280,
        bargap: 0.55,
        shapes: shapes,
        annotations: annotations,
        margin: { t: 26, r: 18, b: 56, l: 56 },
        xaxis: axis({ title: "" }),
        yaxis: axis({
          title: "effect fraction (0–1)",
          range: [0, Math.max(0.1, ...values.filter((v) => Number.isFinite(v)), expected || 0) * 1.2],
        }),
      }
    );
    const legend = $("#legend-mixture");
    if (legend) {
      legend.innerHTML =
        `<span class="item">bars = simulated effect fraction</span>` +
        `<span class="item">hairline = the ${esc(syn.model || "bliss")} expectation` +
        `${expected === null ? "" : " (" + num(expected, 3) + ")"} for two independent agents</span>` +
        `<span class="item">${chip("PREDICTION", "claim:readout")}</span>`;
    }
    if (iso && (iso.points || []).length) {
      /* Normalised isobologram: each axis is the compound's own iso-effect dose,
         so the additivity line is the unit diagonal whatever the two potencies
         are. Raw molar axes put an 8e-11 M and a 3e-8 M compound on scales three
         decades apart and the line reads as flat. */
      const da = Number(iso.d_a) || 1;
      const db = Number(iso.d_b) || 1;
      const pts = iso.points;
      draw(
        "plot-isobologram",
        [
          {
            type: "scatter",
            mode: "lines",
            name: "Loewe additivity",
            x: [1, 0],
            y: [0, 1],
            line: { color: c.base, width: 1 },
            hoverinfo: "skip",
          },
          {
            type: "scatter",
            mode: "markers",
            name: "iso-effect pairs",
            x: pts.map((p) => p.conc_a_M / da),
            y: pts.map((p) => p.conc_b_M / db),
            customdata: pts.map((p) => [p.conc_a_M, p.conc_b_M, p.combination_index]),
            marker: { color: cssVar("--seq-550"), size: 11, line: { width: 2, color: c.surface } },
            hovertemplate:
              "A %{customdata[0]:.2e} M<br>B %{customdata[1]:.2e} M<br>CI %{customdata[2]:.2f}<extra></extra>",
          },
        ],
        {
          height: 280,
          showlegend: true,
          margin: { t: 34, r: 18, b: 52, l: 64 },
          xaxis: axis({ title: `${esc(data.a)} / its own dose`, range: [-0.05, 1.3] }),
          yaxis: axis({ title: `${esc(data.b)} / its own dose`, range: [-0.05, 1.3] }),
        }
      );
      const il = $("#legend-isobologram");
      if (il) {
        il.innerHTML =
          `<span class="item">effect fraction ${num(iso.effect_frac, 2)} on ${esc(iso.readout || iso.assay)}</span>` +
          `<span class="item">A alone = ${sci(iso.d_a, 2)} M · B alone = ${sci(iso.d_b, 2)} M</span>` +
          `<span class="item">combination index = 1 on the line</span>` +
          `<span class="item">${chip("MODEL-DERIVED", "claim:engagement_transformation")}</span>`;
      }
    } else {
      const node = $("#plot-isobologram");
      if (node) node.innerHTML = '<div class="empty">No iso-effect pair exists for this combination.</div>';
    }
    table(
      "tbl-mixture",
      [
        { label: "receptor" },
        { label: "engagement", num: true },
        { label: "model" },
        { label: "components" },
        { label: "not modelled for" },
      ],
      ((mix.mixture || {}).receptors || []).map((row) => [
        esc(row.receptor),
        engagement(row.engagement === undefined ? row.occupancy : row.engagement, 3),
        esc(row.model || ""),
        esc((row.components || []).map((x) => x.compound).join(", ")),
        esc((row.not_modelled_components || []).join(", ") || "—"),
      ]),
      { empty: "No shared receptor." }
    );
  }

  /* mixtures.py states the expectation inside its `why` sentence; parse it so
     the chart's reference line is the library's own number, not a recomputation. */
  function expectedFromWhy(why) {
    const m = /expectation\s+(-?\d+(?:\.\d+)?)/.exec(String(why || ""));
    return m ? Number(m[1]) : null;
  }

  // ==================================================================
  // 0e. EVIDENCE
  // ==================================================================
  const ANALYSES = [
    {
      id: "ablation",
      title: "Model ablation ladder",
      note: "Receptor only, composition only, topology only, full FlyLab — and what each level adds.",
      route: "/api/ablation",
      body: () => ({ compound: design().compound, conc_M: design().conc_M }),
    },
    {
      id: "dependence-landscape",
      title: "Dependence landscape",
      note: "Every compound at every concentration. Minutes of compute; the browser default is the estimate only.",
      route: "/api/dependence/landscape",
      body: () => ({ n: BROWSER_DEPENDENCE_N, estimate_only: false }),
    },
    {
      id: "stability",
      title: "Conclusion stability",
      note: "Every pre-registered conclusion re-derived under every admissible mechanism rule.",
      route: "/api/robustness/stability",
      body: () => ({ fast: true, conc_M: design().conc_M }),
    },
    {
      id: "thresholds",
      title: "Threshold sensitivity",
      note: "The amplify / buffer split recomputed across the whole threshold grid.",
      route: "/api/robustness/thresholds",
      body: () => ({}),
    },
    {
      id: "uncertainty",
      title: "Global uncertainty (Sobol)",
      note: "Variance attribution over the nine uncertain factors. Browser default n_base = 32.",
      route: "/api/uncertainty/global",
      body: () => ({ compound: design().compound, conc_M: design().conc_M, n_base: 32 }),
    },
    {
      id: "voi",
      title: "Value of information",
      note: "Which experiment would remove the most model variance.",
      route: "/api/voi",
      body: () => ({ compound: design().compound, conc_M: design().conc_M, n_base: 32 }),
    },
  ];

  function renderAnalysisCards() {
    const host = $("#analysis-grid");
    if (!host) return;
    host.innerHTML = ANALYSES.map(
      (a) =>
        `<div class="analysis" id="an-${esc(a.id)}"><h4>${esc(a.title)}</h4>` +
        `<p>${esc(a.note)}</p>` +
        `<div class="est" id="est-${esc(a.id)}">estimating…</div>` +
        `<button class="small" type="button" data-analysis="${esc(a.id)}">Run</button></div>`
    ).join("");
    ANALYSES.forEach(async (a) => {
      const node = document.getElementById("est-" + a.id);
      if (!node) return;
      try {
        const est = await post(a.route, Object.assign({}, a.body(), { estimate_only: true }));
        const e = est.runtime_estimate || {};
        node.textContent =
          e.estimate_s === undefined ? "runtime unknown" : `estimated ${e.estimate_s} s before it runs`;
        node.title = e.note || "";
      } catch (err) {
        node.textContent = "estimate unavailable";
      }
    });
  }

  async function runAnalysis(id) {
    const spec = ANALYSES.find((a) => a.id === id);
    const host = $("#analysis-out");
    if (!spec || !host) return;
    host.innerHTML = `<div class="notice info">Running ${esc(spec.title)}…</div>`;
    try {
      const out = await post(spec.route, spec.body());
      host.innerHTML =
        `<h4>${esc(spec.title)}</h4>` +
        `<div class="hint">estimated ${esc((out.runtime_estimate || {}).estimate_s)} s · actual ` +
        `${esc(out.runtime_s === undefined ? "—" : Number(out.runtime_s).toFixed(1))} s</div>` +
        `<pre class="json">${esc(JSON.stringify(out, null, 2).slice(0, 60000))}</pre>`;
      (out.warnings || []).slice(0, 4).forEach((w) => toast(w, "info"));
    } catch (err) {
      host.innerHTML = `<div class="notice">${esc(spec.title)}: ${esc(err.message)}</div>`;
    }
  }

  async function runEvidence() {
    if (!state.cache.dashboard) {
      await runDashboard();
      state.loaded.dashboard = true;
    }
    let validation = state.cache.validation;
    if (!validation) {
      try {
        validation = await api("/api/validation");
        state.cache.validation = validation;
      } catch (err) {
        validation = null;
      }
    }
    renderEvidence();
    return null;
  }

  function renderEvidence() {
    const out = state.cache.dashboard;
    if (!out) return;
    indexProvenance(out);
    const c = colors();
    const key = $("#legend-chips");
    if (key) {
      key.innerHTML = CHIPS.map(
        (k) => `<span class="item">${chip(k)} ${esc(CHIP_NOTE[k])}</span>`
      ).join("");
    }
    table(
      "tbl-evidence",
      [
        { label: "receptor" },
        { label: "engagement", num: true },
        { label: "param" },
        { label: "value (M)", num: true },
        { label: "relation" },
        { label: "species / preparation" },
        { label: "classification" },
      ],
      (out.evidence || []).map((r) => [
        `<span class="swatch" style="background:${r.organism === "insect" ? c.insect : c.vertebrate}"></span>${esc(r.receptor)}`,
        engagement(r.engagement, 3),
        esc(r.param_type || ""),
        sci(r.param_value_M, 2),
        esc(r.relation || ""),
        esc(r.species || ""),
        chip(r.classification, "receptor:" + r.receptor),
      ]),
      { empty: "Run the dashboard first." }
    );

    const claims = out.claims || {};
    const fiu = claims.fact_inference_unknown || {};
    const host = $("#fiu");
    if (host) {
      const col = (title, items, kind) =>
        `<section><h4>${esc(title)} ${chip(kind)}</h4><ul>` +
        (items || []).map((i) => `<li>${esc(i.statement)}</li>`).join("") +
        (items && items.length ? "" : "<li>—</li>") +
        "</ul></section>";
      host.innerHTML =
        col("Fact", fiu.facts, "LITERATURE") +
        col("Model inference", fiu.model_inference, "PREDICTION") +
        col("Unknown", fiu.unknown, "NOT MODELLED");
    }
    const chain = $("#claim-chain");
    if (chain) {
      chain.innerHTML =
        '<div class="chain">' +
        (claims.chain || [])
          .map(
            (l) =>
              `<div class="link"><span class="n">${l.order}</span>` +
              `<span class="step">${esc(l.step)}<br>${chip(l.classification, "claim:" + l.step)}</span>` +
              `<span class="what">${esc(l.statement)}</span></div>`
          )
          .join("") +
        "</div>" +
        `<div class="hint">Labels used by the machine-readable audit: ${esc((claims.labels || []).join(", "))}. ` +
        `Counts for this run: ${esc(JSON.stringify(claims.label_counts || {}))}.</div>`;
    }
    const lib = (out.coverage || {}).library || {};
    const rows = [];
    Object.keys(lib.by_param_type || {}).forEach((k) =>
      rows.push(["parameter type", esc(k), lib.by_param_type[k]])
    );
    Object.keys(lib.by_evidence_tier || {}).forEach((k) =>
      rows.push(["evidence tier", esc(k.replace(/_/g, " ")), lib.by_evidence_tier[k]])
    );
    Object.keys(lib.by_engagement_model || {}).forEach((k) =>
      rows.push(["engagement model", esc(k.replace(/_/g, " ")), lib.by_engagement_model[k]])
    );
    table(
      "tbl-library",
      [{ label: "dimension" }, { label: "value" }, { label: "rows", num: true }],
      rows,
      { empty: "No library census." }
    );

    const val = state.cache.validation || {};
    const vrows = [];
    Object.keys(val.assays || {}).forEach((assayName) => {
      ((val.assays[assayName] || {}).rows || []).forEach((r) =>
        vrows.push(Object.assign({ model_assay: assayName }, r))
      );
    });
    table(
      "tbl-validation",
      [
        { label: "published dataset" },
        { label: "model assay" },
        { label: "compounds", num: true },
        { label: "spearman", num: true },
        { label: "status" },
      ],
      vrows.map((r) => [
        esc(r.id || "—"),
        esc(r.model_assay || "—"),
        esc(r.n_compounds === undefined ? "—" : r.n_compounds),
        r.skipped ? "—" : num(r.spearman_rho, 3),
        r.skipped
          ? `<span class="badge">skipped</span>`
          : r.exact_match
          ? '<span class="badge tier-literature_order">exact order</span>'
          : '<span class="badge tier-class_placeholder">partial</span>',
      ]),
      { empty: "Rank validation is not available in this build." }
    );
    const vhost = $("#tbl-validation").closest(".card");
    let vnote = vhost && vhost.querySelector(".validation-note");
    if (vhost && !vnote) {
      vnote = document.createElement("div");
      vnote.className = "hint validation-note";
      vhost.appendChild(vnote);
    }
    if (vnote) {
      const sum = val.summary || {};
      const kd = val.known_discrepancies || {};
      const n = (Array.isArray(kd) ? kd : kd.discrepancies || []).length;
      vnote.textContent =
        `${sum.n_evaluated || 0} of ${sum.n_runs || 0} runs were scorable (mean rho ` +
        `${sum.mean_rho === undefined ? "—" : Number(sum.mean_rho).toFixed(3)}); ` +
        `${n} literature-vs-library contradictions are on record and are reported, not silently fixed.`;
    }
    renderAnalysisCards();
  }

  // ==================================================================
  // tabs
  // ==================================================================
  const TABS = {
    dashboard: { run: runDashboard, render: renderDashboard },
    dose: { run: runDose, render: renderDose },
    genotype: { run: runGenotype, render: renderGenotype },
    mixtures: { run: runMixtures, render: renderMixtures },
    evidence: { run: runEvidence, render: renderEvidence },
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

  /* Expensive analyses and user-authored experiment designs never start on tab selection. */
  const MANUAL_TABS = ["controls", "experiment"];

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
    const mobilePicker = $("#mobile-panel-picker");
    if (mobilePicker) mobilePicker.value = name;
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
    const validation = state.cache.validation; // library-wide, dose-independent
    state.cache = validation ? { validation: validation } : {};
    state.loaded = {};
    state.dirty = {};
    state.prov = {};
  }

  // ==================================================================
  // wiring
  // ==================================================================
  function wire() {
    $$('[role="tab"]').forEach((b) =>
      b.addEventListener("click", () => activate(b.id.replace("tab-", "")))
    );
    const mobilePicker = $("#mobile-panel-picker");
    if (mobilePicker) mobilePicker.addEventListener("change", () => activate(mobilePicker.value));
    $$("button[data-run]").forEach((b) =>
      b.addEventListener("click", () => runTab(b.getAttribute("data-run")))
    );
    $("#btn-run").addEventListener("click", () => runTab(state.activeTab));

    $("#btn-preset-export").addEventListener("click", () => {
      try {
        const preset = presetForExport();
        download("flylab_design_preset_v1.json", JSON.stringify(preset, null, 2), "application/json;charset=utf-8");
        presetStatus("Exported FlyLab design preset version 1.");
      } catch (err) {
        presetStatus("Could not export: " + err.message);
        toast("preset: " + err.message, "error");
      }
    });
    $("#btn-preset-import").addEventListener("click", () => $("#f-preset-file").click());
    $("#f-preset-file").addEventListener("change", async (event) => {
      const file = event.target.files && event.target.files[0];
      if (!file) return;
      try {
        if (file.size > 1024 * 1024) throw new Error("preset file is larger than 1 MB");
        const parsed = JSON.parse(await file.text());
        applyPreset(parsed);
        presetStatus("Imported preset. Open Experiment and choose Run design to calculate results.");
        toast("design preset imported", "ok");
      } catch (err) {
        const message = (err && err.message) || "could not read preset file";
        presetStatus("Import failed: " + message);
        toast("preset: " + message, "error");
      } finally {
        event.target.value = "";
      }
    });

    $("#f-logconc").addEventListener("input", (e) => {
      setConc(Math.pow(10, Number(e.target.value)));
      invalidate();
    });
    $("#f-logconc").addEventListener("change", () => persistDesign());
    $("#f-conc").addEventListener("change", (e) => {
      const v = Number(e.target.value);
      if (Number.isFinite(v) && v > 0) {
        setConc(v);
        persistDesign();
      } else {
        toast("concentration must be a positive number in molar", "error");
        setConc(Math.pow(10, Number($("#f-logconc").value)));
      }
      invalidate();
    });
    $("#f-compound").addEventListener("change", () => {
      updateCompoundBadges();
      invalidate();
      persistDesign();
    });
    $("#f-graph").addEventListener("change", () => {
      updateGraphHint();
      $("#f-cy-minw").value = $("#f-graph").value === "taste_motor" ? 12 : 5;
      invalidate();
      persistDesign();
    });
    ["#f-engine", "#f-seed", "#f-sugar", "#f-bitter", "#f-genotype"].forEach((sel) =>
      $(sel).addEventListener("change", () => {
        invalidate();
        persistDesign();
      })
    );

    /* The dose rail drives every dashboard card. `input` only moves the label
       (so dragging stays smooth); the run happens on `change`, i.e. release. */
    const rail = $("#f-dash-logconc");
    if (rail) {
      rail.addEventListener("input", (e) => {
        const v = Math.pow(10, Number(e.target.value));
        const read = $("#dash-conc-read");
        if (read) read.textContent = sci(v) + " M";
      });
      rail.addEventListener("change", (e) => {
        setConc(Math.pow(10, Number(e.target.value)));
        invalidate();
        persistDesign();
        runTab(state.activeTab);
      });
    }
    $$("#dash-ladder button").forEach((b) =>
      b.addEventListener("click", () => {
        setConc(Number(b.getAttribute("data-conc")));
        invalidate();
        persistDesign();
        runTab(state.activeTab);
      })
    );

    const cmp = $("#f-cmp-compounds");
    if (cmp) cmp.addEventListener("change", () => estimateCompare());
    const cmpDep = $("#f-cmp-dependence");
    if (cmpDep) cmpDep.addEventListener("change", () => estimateCompare());
    const cmpRun = $("#btn-cmp-run");
    if (cmpRun) {
      cmpRun.addEventListener("click", () => {
        cmpRun.disabled = true;
        status("busy", "running compare…");
        runCompare()
          .then((ms) => status("ok", "compare ready", ms))
          .catch((err) => {
            status("error", "compare failed");
            toast("compare: " + err.message, "error");
          })
          .finally(() => {
            cmpRun.disabled = false;
          });
      });
    }

    const genoOther = $("#btn-geno-other");
    if (genoOther) {
      genoOther.addEventListener("click", () =>
        runGenotypeOther().catch((err) => toast("genotype: " + err.message, "error"))
      );
    }

    const provClose = $("#btn-prov-close");
    if (provClose) provClose.addEventListener("click", closeProvenance);
    document.addEventListener("keydown", (ev) => {
      if (ev.key === "Escape") {
        closeProvenance();
        tourEnd();
      }
      if (!tour.open) return;
      if (ev.key === "ArrowRight") tourGo(tour.i + 1);
      if (ev.key === "ArrowLeft") tourGo(tour.i - 1);
    });

    ["#btn-tour", "#btn-tour-inline"].forEach((sel) => {
      const btn = $(sel);
      if (btn) btn.addEventListener("click", () => (tour.open ? tourEnd() : tourStart(true)));
    });
    window.addEventListener("scroll", tourSchedulePlace, { passive: true });
    window.addEventListener("resize", tourSchedulePlace);

    $("#btn-ic50").addEventListener("click", () =>
      runIc50().catch((err) => toast("IC50: " + err.message, "error"))
    );
    $("#f-cy-state").addEventListener("change", renderCircuit);
    $("#btn-cy-fit").addEventListener("click", () => {
      if (state.cy) state.cy.fit(undefined, 24);
    });

    document.addEventListener("click", (ev) => {
      const chipEl = ev.target.closest && ev.target.closest(".chip[data-prov]");
      if (chipEl) {
        openProvenance(chipEl.getAttribute("data-prov"));
        return;
      }
      const chipWhat = ev.target.closest && ev.target.closest(".chip[data-explain]");
      if (chipWhat) {
        explainChip(chipWhat.getAttribute("data-explain"));
        return;
      }
      const tourBtn = ev.target.closest && ev.target.closest("[data-tour]");
      if (tourBtn) {
        const act = tourBtn.getAttribute("data-tour");
        if (act === "next") tourGo(tour.i + 1);
        else if (act === "prev") tourGo(tour.i - 1);
        else tourEnd();
        return;
      }
      const an = ev.target.closest && ev.target.closest("[data-analysis]");
      if (an) {
        an.disabled = true;
        runAnalysis(an.getAttribute("data-analysis")).finally(() => {
          an.disabled = false;
        });
        return;
      }
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

    ["#f-exp-start", "#f-exp-stop", "#f-exp-ppd"].forEach((sel) => {
      $(sel).addEventListener("input", updateLadderHint);
      $(sel).addEventListener("change", () => persistDesign({ experimentSettings: experimentSettings() }));
    });
    ["#f-exp-compounds", "#f-exp-readouts", "#f-exp-assay", "#f-exp-reps", "#f-exp-vehicle"].forEach((sel) => {
      $(sel).addEventListener("change", () => persistDesign({ experimentSettings: experimentSettings() }));
    });
    $("#f-exp-blind").addEventListener("change", () => {
      state.revealed = false;
      $("#btn-exp-reveal").disabled = !$("#f-exp-blind").checked || !state.experiment;
      renderExperiment();
      persistDesign({ experimentSettings: experimentSettings() });
    });
    $("#f-exp-randomize").addEventListener("change", () => {
      state.runOrder = null;
      renderExperiment();
      persistDesign({ experimentSettings: experimentSettings() });
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
  // guided tour
  // ==================================================================
  /* A ring around a real panel and a card that explains it. Two rules shape
     the whole thing:

       * it never blocks the bench. There is no backdrop and the ring carries
         pointer-events:none, so every control stays clickable while the tour
         is open. A tour that has to be dismissed before the tool can be used
         is a dialog, not a tour.
       * every step states one interpretation caveat, because the failure mode
         this bench has is not "I cannot find the button", it is "I read the
         number as something it is not".

     No library, no build step, no network: the same file serves the FastAPI
     bench and the Pyodide build. Dismissal is remembered per viewer in
     localStorage (see store(), which survives blocked storage). */
  const TOUR_VERSION = 1;
  const TOUR_GUIDE_URL =
    "https://github.com/pinkysworld/FlyLab/blob/main/docs/INTERPRETATION.md";

  const TOUR_STEPS = [
    {
      sel: "#dash-overview",
      title: "Compound, dose, and what is missing",
      what: "The header names the compound, the dose everything below is computed at, and how many receptor rows carry a sourced value.",
      caveat:
        "The rows counted as “not modelled” have no number at all and are excluded from every figure on this page — they are not zeros.",
    },
    {
      sel: "#dash-tiles",
      title: "The four headline tiles",
      what: "Insect engagement, vertebrate engagement, the potency ratio between them, and the evidence tier the whole row rests on.",
      caveat:
        "Engagement is not occupancy: a tile can read 1.000 from a binding constant measured in another species. Click a chip for the parameter type and the evidence distance.",
    },
    {
      sel: ".dose-rail",
      title: "Move the dose",
      what: "Every card recomputes at the concentration you release the slider on; the ladder buttons jump to round decades.",
      caveat:
        "Above saturation the curve flattens because the receptor is full, not because the model is confident. At 1 µM a two-fold error in imidacloprid's potency moves the circuit by 0.000 Hz.",
    },
    {
      sel: "#plot-dash-circuit",
      title: "Circuit consequence",
      what: "Vehicle against treated firing rates for named cells on the committed MaleCNS cut.",
      caveat:
        "A relative change in a simulated rate. There is no dose, no exposure and no living fly anywhere behind these numbers.",
    },
    {
      sel: "#dash-dep-verdict",
      title: "Is the effect wiring-dependent?",
      what: "Permutation nulls destroy one kind of network structure at a time, and the verdict names the weakest one the test cannot tell apart from the real cut.",
      caveat:
        "A large permutation probability licenses only “not distinguishable at n shuffles”. And the verdict belongs to the cut, and to its size: imidacloprid is composition-dominated on the 1126-cell `named` in-star and topology-dependent on `taste_motor`, and it moves again on larger cuts of the same connectome. Never quote a verdict without naming the cut and its size.",
    },
    {
      sel: "#card-trust",
      title: "What can I trust?",
      what: "Each layer of the chain reported separately: sourced values, asserted rules, measured wiring, predicted transmitters, and the validation that does not exist.",
      caveat:
        "There is deliberately no blended confidence score. The layers fail independently, and one number would hide which of them is weak.",
    },
    {
      sel: "#card-compare",
      title: "Compare compounds",
      what: "Receptor selectivity beside circuit selectivity for a set of compounds at the same dose.",
      caveat:
        "The two indices rank compounds differently, and neither is a safety margin — a ratio on a teaching library says nothing about a vertebrate.",
    },
  ];

  const tour = { open: false, i: 0, ring: null, card: null, raf: null, restore: null };

  function tourEls() {
    if (!tour.ring) {
      tour.ring = document.createElement("div");
      tour.ring.className = "tour-ring";
      tour.ring.setAttribute("aria-hidden", "true");
      document.body.appendChild(tour.ring);
    }
    if (!tour.card) {
      tour.card = document.createElement("div");
      tour.card.className = "tour-card";
      tour.card.id = "tour-card";
      tour.card.setAttribute("role", "dialog");
      // aria-modal stays false on purpose: the page behind is still live
      tour.card.setAttribute("aria-modal", "false");
      tour.card.setAttribute("aria-labelledby", "tour-title");
      document.body.appendChild(tour.card);
    }
    return tour;
  }

  function tourTarget(i) {
    const step = TOUR_STEPS[i];
    if (!step) return null;
    const el = $(step.sel);
    if (!el) return null;
    // a card that has not been run yet is still worth pointing at; an element
    // with no box at all (display:none) is not
    const box = el.getBoundingClientRect();
    if (!box.width && !box.height) return null;
    return el.closest(".card") || el;
  }

  function tourPlace() {
    if (!tour.open) return;
    const el = tourTarget(tour.i);
    const ring = tour.ring;
    const card = tour.card;
    if (!ring || !card) return;
    const sx = window.scrollX || window.pageXOffset || 0;
    const sy = window.scrollY || window.pageYOffset || 0;
    if (!el) {
      // the panel is not on screen (another tab, or a build without it):
      // centre the card and hide the ring rather than dropping the step
      ring.style.display = "none";
      card.style.top = sy + Math.max(16, window.innerHeight / 2 - 120) + "px";
      card.style.left = sx + Math.max(12, (window.innerWidth - card.offsetWidth) / 2) + "px";
      return;
    }
    const b = el.getBoundingClientRect();
    ring.style.display = "block";
    ring.style.top = b.top + sy - 4 + "px";
    ring.style.left = b.left + sx - 4 + "px";
    ring.style.width = b.width + 8 + "px";
    ring.style.height = b.height + 8 + "px";

    const cw = card.offsetWidth || 360;
    const ch = card.offsetHeight || 180;
    // below the target when there is room, otherwise above it, and never off
    // the left or right edge of the viewport
    let top = b.bottom + sy + 12;
    if (b.bottom + ch + 24 > window.innerHeight && b.top - ch - 12 > 0) {
      top = b.top + sy - ch - 12;
    }
    let left = b.left + sx;
    const maxLeft = sx + window.innerWidth - cw - 14;
    if (left > maxLeft) left = maxLeft;
    if (left < sx + 12) left = sx + 12;
    card.style.top = Math.max(sy + 8, top) + "px";
    card.style.left = left + "px";
  }

  function tourSchedulePlace() {
    if (tour.raf) cancelAnimationFrame(tour.raf);
    tour.raf = requestAnimationFrame(() => {
      tour.raf = null;
      tourPlace();
    });
  }

  function tourRender() {
    const step = TOUR_STEPS[tour.i];
    if (!step) return;
    const dots = TOUR_STEPS.map(
      (_, k) => `<i data-now="${k === tour.i ? "true" : "false"}"></i>`
    ).join("");
    const last = tour.i === TOUR_STEPS.length - 1;
    tour.card.innerHTML =
      `<p class="tour-step">Step ${tour.i + 1} of ${TOUR_STEPS.length}</p>` +
      `<h3 id="tour-title">${esc(step.title)}</h3>` +
      `<p>${esc(step.what)}</p>` +
      `<p class="tour-caveat">${esc(step.caveat)}</p>` +
      `<div class="tour-actions">` +
      `<div class="tour-dots" aria-hidden="true">${dots}</div>` +
      `<span class="spacer"></span>` +
      (tour.i > 0 ? `<button type="button" class="small" data-tour="prev">Back</button>` : "") +
      `<button type="button" class="small" data-tour="end">Dismiss</button>` +
      `<button type="button" class="primary small" data-tour="${last ? "end" : "next"}">` +
      (last ? "Done" : "Next") +
      `</button>` +
      `</div>` +
      `<p style="margin:8px 0 0"><a class="howto" href="${TOUR_GUIDE_URL}" target="_blank" rel="noreferrer noopener">How to read this &rarr;</a></p>`;
    const el = tourTarget(tour.i);
    if (el && el.scrollIntoView) {
      try {
        el.scrollIntoView({ block: "center", behavior: "smooth" });
      } catch (err) {
        el.scrollIntoView();
      }
    }
    tourSchedulePlace();
    // one reflow later the smooth scroll has landed
    setTimeout(tourSchedulePlace, 420);
  }

  function tourGo(i) {
    if (i < 0 || i >= TOUR_STEPS.length) return tourEnd();
    tour.i = i;
    tourRender();
  }

  function tourStart(force) {
    if (!force && store("tourSeen") === TOUR_VERSION) return;
    if (state.activeTab !== "dashboard") activate("dashboard");
    tourEls();
    tour.open = true;
    tour.i = 0;
    tour.restore = document.activeElement;
    tour.card.hidden = false;
    tour.ring.hidden = false;
    const btn = $("#btn-tour");
    if (btn) btn.setAttribute("aria-expanded", "true");
    tourRender();
    // focus the card so keyboard users land on the controls, but only after
    // the scroll has settled; the page keeps working either way
    setTimeout(() => {
      const first = tour.card && tour.card.querySelector("button");
      if (tour.open && first) first.focus({ preventScroll: true });
    }, 440);
  }

  function tourEnd() {
    if (!tour.open) return;
    tour.open = false;
    store("tourSeen", TOUR_VERSION);
    if (tour.ring) tour.ring.hidden = true;
    if (tour.card) tour.card.hidden = true;
    const btn = $("#btn-tour");
    if (btn) btn.setAttribute("aria-expanded", "false");
    const back = tour.restore;
    tour.restore = null;
    if (back && back.focus) {
      try {
        back.focus({ preventScroll: true });
      } catch (err) {
        /* the element may be gone after a re-render */
      }
    }
  }

  /* A classification chip with no provenance record still has something to
     say: what the label itself means. It opens the same drawer, so there is
     one place a reader looks for "where did this come from". */
  function explainChip(kind) {
    const k = CHIPS.indexOf(kind) >= 0 ? kind : "NOT MODELLED";
    const drawer = $("#provenance");
    const body = $("#prov-body");
    const title = $("#prov-title");
    if (!drawer || !body) return;
    if (title) title.textContent = "What this label means";
    body.innerHTML =
      `<div class="prov-head"><span class="name">${esc(k)}</span>${chip(k, null, k.toLowerCase())}</div>` +
      `<p>${esc(CHIP_NOTE[k])}</p>` +
      `<dl><dt>the whole vocabulary</dt><dd><ul>` +
      CHIPS.map((c) => `<li><b>${esc(c)}</b> — ${esc(CHIP_NOTE[c])}</li>`).join("") +
      `</ul></dd></dl>` +
      `<p>Chips carrying a provenance record open the source instead. ` +
      `<a href="${TOUR_GUIDE_URL}" target="_blank" rel="noreferrer noopener">How to read this &rarr;</a></p>`;
    drawer.dataset.open = "true";
    drawer.setAttribute("aria-hidden", "false");
  }

  // ==================================================================
  // static build
  // ==================================================================
  /* flylab.circuit.lif.BROWSER_T_MS: the browser keeps the 0.1 ms integration
     step and shortens the window instead, so a run here is the same model as
     a run on the server, not a coarser one. */
  const STATIC_T_MS = 200;

  function staticNote(panelId, text) {
    const head = document.querySelector("#panel-" + panelId + " .panel-head");
    if (!head || head.querySelector(".static-note")) return;
    const p = document.createElement("p");
    p.className = "sub static-note";
    p.textContent = text;
    const heading = head.querySelector("h1");
    if (heading && heading.nextSibling) head.insertBefore(p, heading.nextSibling);
    else head.appendChild(p);
  }

  function applyStaticDefaults() {
    const tms = $("#f-tms");
    if (tms) tms.value = String(STATIC_T_MS);
    staticNote(
      "spikes",
      `Browser build: ${STATIC_T_MS} ms window by default (same 0.1 ms step as the served bench, shorter run), ` +
        "and the rate engine is the default elsewhere."
    );
    const build = window.FLYLAB_BUILD || {};
    toast(
      `static build: Python ${build.pyodide_version ? "(Pyodide " + build.pyodide_version + ") " : ""}runs in this tab; nothing is uploaded`,
      "info"
    );
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
    if (isStatic()) {
      try {
        await window.flylabReady;
      } catch (err) {
        status("error", "the in-browser Python runtime did not start");
        return;
      }
      applyStaticDefaults();
    }
    if (!hasPlotly()) toast("Plotly did not load from either CDN; charts fall back to tables", "info");
    if (!hasCytoscape()) toast("cytoscape.js did not load from either CDN; the circuit view falls back to tables", "info");

    try {
      await loadMeta();
    } catch (err) {
      status("error", "could not reach the bench server");
      toast("meta: " + err.message, "error");
      return;
    }

    const saved = store("design") || {};
    restoreExperimentSettings(saved);
    $("#f-exp-reps").value = $("#f-exp-reps").value || 2;

    estimateCompare().catch(() => {});
    status("ok", "ready");
    activate("dashboard");
    /* first visit only: the tour remembers its own dismissal, and it never
       blocks the bench, so a returning viewer never sees it again unless the
       "Guided tour" button asks for it. */
    setTimeout(() => {
      try {
        tourStart(false);
      } catch (err) {
        /* a tour must never be the reason the bench fails to load */
      }
    }, 900);
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", boot);
  else boot();
})();
