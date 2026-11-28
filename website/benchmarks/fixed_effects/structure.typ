// Figures of the data generating processes of the fixed-effects benchmark
// (see ../fixed_effects.typ), one per process: `dgp-figure("balanced")`,
// `"random"`, `"round-robin"` and `"low-mobility"`. Each is a small grid:
// the years every individual (row) spends in each year or group (column),
// labelled with an example of its process (states and years, shoppers and
// products, students and teachers, workers and firms).
// The data are simulated here (a small linear congruential generator), at
// a reduced size (see `dgp-canvas`), over 10 years.
#import "@preview/cetz:0.5.2"

#let _n-year = 10

#let _next(s) = calc.rem(s * 1103515245 + 12345, 2147483648)
#let _int(s, m) = calc.rem(calc.quo(s, 65536), m)
#let _zeros(nr, nc) = range(nr).map(_ => range(nc).map(_ => 0))

// The group (firm, product, teacher) of each of `ni` individuals in each
// year: one array of `_n-year` groups, out of `nf`, per individual
#let _groups(type, ni, nf) = {
  let s = 42
  let all = ()
  for i in range(ni) {
    // low mobility: individuals start at a firm (sorted, for display)
    let f = calc.quo(i * nf, ni)
    let own = ()
    for t in range(_n-year) {
      if type == "random" {
        s = _next(s)
        f = _int(s, nf)
      } else if type == "round-robin" {
        // round robin over the rows, as in dgp.R
        f = calc.rem(i * _n-year + t, nf)
      } else {
        // they move to a random firm with probability 0.08 a year
        s = _next(s)
        if t > 0 and _int(s, 100) < 8 {
          s = _next(s)
          f = _int(s, nf)
        }
      }
      own.push(f)
    }
    all.push(own)
  }
  all
}

// individuals x groups: the years of each individual in each group
#let _crosstab(groups, nf) = groups.map(own => range(nf).map(f => own.filter(g => g == f).len()))

#let _ramp = gradient.linear(rgb("#a9d5cd"), rgb("#0b4a45"))
#let _dark = rgb("#0b4a45")
#let _empty = luma(128, 14%)
#let _ink = luma(115)
#let _label(body) = text(8.5pt, fill: _ink, body)

// One heatmap with its top left corner at (x, y), `w` by `h`: shades from
// 1 (light) to the number of years (dark), zeros left empty.
#let _heat(x, y, w, h, M) = {
  import cetz.draw: *
  let nr = M.len()
  let nc = M.at(0).len()
  let cw = w / nc
  let ch = h / nr
  rect((x, y), (x + w, y - h), fill: _empty, stroke: none)
  for r in range(nr) {
    for j in range(nc) {
      let v = M.at(r).at(j)
      if v > 0 {
        rect((x + j * cw + 0.04 * cw, y - r * ch - 0.04 * ch),
             (x + (j + 1) * cw - 0.04 * cw, y - (r + 1) * ch + 0.04 * ch),
             fill: _ramp.sample((v - 1) * 100% / (_n-year - 1)), stroke: none)
      }
    }
  }
}

// The shades of the heatmap, from 1 to 10 years
#let _legend(x, y, what) = {
  import cetz.draw: *
  let x0 = x + 3.3
  content((x0 - 0.12, y), _label[#what: 1], anchor: "east")
  for k in range(10) {
    rect((x0 + k * 0.22, y - 0.09), (x0 + (k + 1) * 0.22 - 0.02, y + 0.09),
         fill: _ramp.sample(k * 100% / 9), stroke: none)
  }
  content((x0 + 2.3, y), _label[10], anchor: "west")
}

#let dgp-canvas(type) = cetz.canvas(length: 1cm, {
  import cetz.draw: *
  let balanced = type == "balanced"
  let chain = type == "round-robin"
  // individuals (rows) and groups (columns) drawn: a small grid
  let (ni, nf) = (
    "balanced": (12, _n-year),
    "random": (12, 24),
    "round-robin": (9, 45),
    "low-mobility": (12, 24),
  ).at(type)
  let groups = if balanced { range(ni).map(_ => range(nf)) } else { _groups(type, ni, nf) }
  // round robin: students in the order of their first teacher (stairs)
  if chain { groups = groups.sorted(key: own => own.at(0)) }
  let N = _crosstab(groups, nf)
  let cw = if chain { 0.2 } else { 0.26 }   // cell width and height
  let ch = 0.3
  let (w, h) = (nf * cw, ni * ch)
  // the example shown: who (rows), what they meet (columns)
  let (who, what) = (
    "balanced": ([states], [years]),
    "random": ([shoppers], [products]),
    "round-robin": ([students], [teachers]),
    "low-mobility": ([workers], [firms]),
  ).at(type)

  _heat(0, 0, w, h, N)
  content((-0.3, -h / 2), angle: 90deg, _label[#ni #who])
  content((w / 2, -h - 0.3), _label[#nf #what])
  if balanced {
    for t in range(nf) {
      content(((t + 0.5) * cw, 0.25), text(7pt, fill: _ink)[#(t + 1)])
    }
  }
  if type == "random" { _legend(0, -h - 0.85, [Purchases of the product]) }
  if type == "low-mobility" { _legend(0, -h - 0.85, [Years at the firm]) }
})

// On the website (HTML) the drawing is an inline SVG frame.
#let dgp-figure(type) = context {
  if target() == "html" { html.frame(dgp-canvas(type)) } else { dgp-canvas(type) }
}
