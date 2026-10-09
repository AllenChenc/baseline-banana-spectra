let rawX = [];
let rawY = [];

let processedX = [];
let processedY = [];
let baselineY = [];

let currentFileName = "spectrum";

const fileInput = document.getElementById("fileInput");
const fileName = document.getElementById("fileName");
const pointCount = document.getElementById("pointCount");
const statusText = document.getElementById("status");

const reverseX = document.getElementById("reverseX");
const baselineMethod = document.getElementById("baselineMethod");
const polyDegree = document.getElementById("polyDegree");
const normalizeMethod = document.getElementById("normalizeMethod");

const applyButton = document.getElementById("applyButton");
const resetButton = document.getElementById("resetButton");
const exportButton = document.getElementById("exportButton");


/* -------------------------
   FILE IMPORT
------------------------- */

fileInput.addEventListener("change", function(event) {

  const file = event.target.files[0];

  if (!file) return;

  currentFileName = file.name.replace(/\.[^/.]+$/, "");

  const reader = new FileReader();

  reader.onload = function(e) {

    const text = e.target.result;

    const lines = text
      .split(/\r?\n/)
      .filter(line => line.trim() !== "");

    const x = [];
    const y = [];

    for (const line of lines) {

      const parts = line
        .trim()
        .split(/[\s,;\t]+/);

      if (parts.length < 2) continue;

      const xv = parseFloat(parts[0]);
      const yv = parseFloat(parts[1]);

      if (!Number.isNaN(xv) && !Number.isNaN(yv)) {
        x.push(xv);
        y.push(yv);
      }
    }

    if (x.length < 2) {
      statusText.textContent =
        "Could not detect valid two-column numeric data.";
      return;
    }

    rawX = x;
    rawY = y;

    processedX = [...rawX];
    processedY = [...rawY];
    baselineY = [];

    fileName.textContent = file.name;
    pointCount.textContent =
      `${rawX.length} data points`;

    statusText.textContent =
      "Raw spectrum loaded.";

    baselineMethod.value = "none";
    normalizeMethod.value = "none";

    drawPlot();
  };

  reader.readAsText(file);
});


/* -------------------------
   PLOT
------------------------- */

function drawPlot() {

  if (rawX.length === 0) return;

  const traces = [];

  traces.push({
    x: rawX,
    y: rawY,
    type: "scatter",
    mode: "lines",
    name: "Raw",
    line: {
      width: 1.5
    }
  });

  if (baselineY.length === rawY.length) {

    traces.push({
      x: rawX,
      y: baselineY,
      type: "scatter",
      mode: "lines",
      name: "Baseline",
      line: {
        dash: "dash",
        width: 2
      }
    });

  }

  if (
    processedY.length === rawY.length &&
    (
      baselineMethod.value !== "none" ||
      normalizeMethod.value !== "none"
    )
  ) {

    traces.push({
      x: processedX,
      y: processedY,
      type: "scatter",
      mode: "lines",
      name: "Processed",
      line: {
        width: 2
      }
    });

  }

  const layout = {

    margin: {
      t: 40,
      r: 30,
      b: 70,
      l: 80
    },

    xaxis: {
      title: "Wavenumber / Raman Shift",
      autorange:
        reverseX.checked
          ? "reversed"
          : true
    },

    yaxis: {
      title: "Intensity"
    },

    legend: {
      orientation: "h",
      x: 0,
      y: 1.08
    },

    hovermode: "closest"
  };

  Plotly.react(
    "plot",
    traces,
    layout,
    {
      responsive: true,
      displaylogo: false
    }
  );
}


reverseX.addEventListener("change", drawPlot);


/* -------------------------
   POLYNOMIAL FIT
------------------------- */

function polynomialFit(x, y, degree) {

  const n = degree + 1;

  const A = Array.from(
    { length: n },
    () => Array(n).fill(0)
  );

  const B = Array(n).fill(0);

  /*
    Scale x before fitting.

    This greatly improves numerical stability
    for spectroscopy axes such as 2800–3800.
  */

  const xmin = Math.min(...x);
  const xmax = Math.max(...x);

  const xs = x.map(v =>
    2 * (v - xmin) / (xmax - xmin) - 1
  );

  for (let row = 0; row < n; row++) {

    for (let col = 0; col < n; col++) {

      let sum = 0;

      for (let i = 0; i < xs.length; i++) {
        sum += Math.pow(xs[i], row + col);
      }

      A[row][col] = sum;
    }

    let sumB = 0;

    for (let i = 0; i < xs.length; i++) {
      sumB += y[i] * Math.pow(xs[i], row);
    }

    B[row] = sumB;
  }

  const coefficients = solveLinearSystem(A, B);

  return xs.map(value => {

    let result = 0;

    for (let j = 0; j < coefficients.length; j++) {
      result +=
        coefficients[j] *
        Math.pow(value, j);
    }

    return result;
  });
}


function solveLinearSystem(A, B) {

  const n = B.length;

  const M = A.map(
    (row, i) => [...row, B[i]]
  );

  for (let i = 0; i < n; i++) {

    /*
      Partial pivoting
    */

    let maxRow = i;

    for (let k = i + 1; k < n; k++) {
      if (
        Math.abs(M[k][i]) >
        Math.abs(M[maxRow][i])
      ) {
        maxRow = k;
      }
    }

    [M[i], M[maxRow]] =
      [M[maxRow], M[i]];

    const pivot = M[i][i];

    if (Math.abs(pivot) < 1e-12) {
      throw new Error(
        "Polynomial fit failed."
      );
    }

    for (let j = i; j <= n; j++) {
      M[i][j] /= pivot;
    }

    for (let k = 0; k < n; k++) {

      if (k === i) continue;

      const factor = M[k][i];

      for (let j = i; j <= n; j++) {
        M[k][j] -= factor * M[i][j];
      }
    }
  }

  return M.map(row => row[n]);
}


/* -------------------------
   NORMALIZATION
------------------------- */

function normalizeData(x, y, method) {

  if (method === "none") {
    return [...y];
  }

  if (method === "max") {

    const maxValue =
      Math.max(...y.map(v => Math.abs(v)));

    if (maxValue === 0) return [...y];

    return y.map(v => v / maxValue);
  }

  if (method === "minmax") {

    const minValue = Math.min(...y);
    const maxValue = Math.max(...y);

    const range = maxValue - minValue;

    if (range === 0) return [...y];

    return y.map(
      v => (v - minValue) / range
    );
  }

  if (method === "area") {

    let area = 0;

    for (let i = 1; i < x.length; i++) {

      const dx =
        Math.abs(x[i] - x[i - 1]);

      area +=
        0.5 *
        (Math.abs(y[i]) +
         Math.abs(y[i - 1])) *
        dx;
    }

    if (area === 0) return [...y];

    return y.map(v => v / area);
  }

  return [...y];
}


/* -------------------------
   APPLY PROCESSING
------------------------- */

applyButton.addEventListener("click", function() {

  if (rawX.length === 0) {

    statusText.textContent =
      "Upload a spectrum first.";

    return;
  }

  try {

    processedX = [...rawX];

    let workingY = [...rawY];

    baselineY = [];

    if (
      baselineMethod.value ===
      "polynomial"
    ) {

      const degree =
        parseInt(polyDegree.value);

      baselineY =
        polynomialFit(
          rawX,
          rawY,
          degree
        );

      workingY =
        rawY.map(
          (value, index) =>
            value - baselineY[index]
        );
    }

    workingY =
      normalizeData(
        processedX,
        workingY,
        normalizeMethod.value
      );

    processedY = workingY;

    statusText.textContent =
      "Processing applied.";

    drawPlot();

  } catch (error) {

    console.error(error);

    statusText.textContent =
      "Processing failed: " +
      error.message;
  }
});


/* -------------------------
   RESET
------------------------- */

resetButton.addEventListener("click", function() {

  if (rawX.length === 0) return;

  processedX = [...rawX];
  processedY = [...rawY];
  baselineY = [];

  baselineMethod.value = "none";
  normalizeMethod.value = "none";
  reverseX.checked = false;

  statusText.textContent =
    "Processing reset.";

  drawPlot();
});


/* -------------------------
   EXPORT CSV
------------------------- */

exportButton.addEventListener("click", function() {

  if (processedX.length === 0) {

    statusText.textContent =
      "Nothing to export.";

    return;
  }

  let csv =
    "X,Processed_Intensity\n";

  for (
    let i = 0;
    i < processedX.length;
    i++
  ) {

    csv +=
      `${processedX[i]},${processedY[i]}\n`;
  }

  const blob =
    new Blob(
      [csv],
      { type: "text/csv;charset=utf-8;" }
    );

  const url =
    URL.createObjectURL(blob);

  const link =
    document.createElement("a");

  link.href = url;

  link.download =
    `${currentFileName}_processed.csv`;

  document.body.appendChild(link);

  link.click();

  document.body.removeChild(link);

  URL.revokeObjectURL(url);
});
