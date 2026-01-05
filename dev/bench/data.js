window.BENCHMARK_DATA = {
  "lastUpdate": 1767607268415,
  "repoUrl": "https://github.com/jipolanco/VortexPasta.jl",
  "entries": {
    "Julia benchmark result": [
      {
        "commit": {
          "author": {
            "email": "juan-ignacio.polanco@cnrs.fr",
            "name": "Juan Ignacio Polanco",
            "username": "jipolanco"
          },
          "committer": {
            "email": "juan-ignacio.polanco@cnrs.fr",
            "name": "Juan Ignacio Polanco",
            "username": "jipolanco"
          },
          "distinct": true,
          "id": "3018e92836841888b340c3b3fa0425cfa04643f6",
          "message": "Minor changes to CI",
          "timestamp": "2026-01-05T10:46:01+01:00",
          "tree_id": "e7ea1fa3afae3f71ae3d9ad67831dc1dadb904f3",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/3018e92836841888b340c3b3fa0425cfa04643f6"
        },
        "date": 1767607261680,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 13144066,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367264\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 13113239,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 793326970,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5317088\nallocs=15937\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 667412525,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5210576\nallocs=15508\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8971907.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=212560\nallocs=2051\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1222057444,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8979632\nallocs=55324\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1530137265.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42098944\nallocs=24956\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5613166,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73488\nallocs=1709\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 880411,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84176\nallocs=2043\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 5347173.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=244912\nallocs=5821\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3166476128,
            "unit": "ns",
            "extra": "gctime=0\nmemory=323521520\nallocs=4179258\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 170746909,
            "unit": "ns",
            "extra": "gctime=17104014\nmemory=169989864\nallocs=5995736\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 7157500,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/set_elements!",
            "value": 2443875,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_pair",
            "value": 95840585,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/iterator_interface",
            "value": 160049439.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=83202144\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_source",
            "value": 117648013,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/set_elements!",
            "value": 9295956,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_pair",
            "value": 357059365.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/iterator_interface",
            "value": 422281147.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_source",
            "value": 340471545,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2989091747,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71446288\nallocs=2241524\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 27820031,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131280\nallocs=2528\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}