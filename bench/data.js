window.BENCHMARK_DATA = {
  "lastUpdate": 1776675729500,
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
          "id": "ed6361e6927c42e7b6ea5b3efa576d3ddc2046c6",
          "message": "Update benchmark.yml",
          "timestamp": "2026-01-05T11:26:28+01:00",
          "tree_id": "43a3c1f2d12e40b0edb4d7670c52c5fcb53db1ff",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/ed6361e6927c42e7b6ea5b3efa576d3ddc2046c6"
        },
        "date": 1767609228274,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 12821542,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367264\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12659828,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 724295609,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5329184\nallocs=16141\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 604192608,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5216624\nallocs=15610\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8638409.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=212560\nallocs=2051\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1162603509,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9034064\nallocs=56242\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1526021727.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42153376\nallocs=25874\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5477070,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73488\nallocs=1709\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 829285,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4016\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 5394225.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=244912\nallocs=5821\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3261084097.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=330872128\nallocs=4199215\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 152859964,
            "unit": "ns",
            "extra": "gctime=16479858\nmemory=169990008\nallocs=5995739\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 6679431,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/set_elements!",
            "value": 2382400,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_pair",
            "value": 61261579.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/iterator_interface",
            "value": 96098658.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=83202144\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_source",
            "value": 61703622,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/set_elements!",
            "value": 7695996,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_pair",
            "value": 300289906,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/iterator_interface",
            "value": 387711864,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_source",
            "value": 299841051,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2811073585,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71476528\nallocs=2242034\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 28045601.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131280\nallocs=2528\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "a00c5c6644e46dd2b842b95346c389a55422f4ff",
          "message": "Try to fix vector_of_vector tests on Julia 1.12",
          "timestamp": "2026-01-05T11:26:41+01:00",
          "tree_id": "5fee92f25a16699ad397d66cc85d43db086daf41",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/a00c5c6644e46dd2b842b95346c389a55422f4ff"
        },
        "date": 1767609248108,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 13050236,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367264\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12852513,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 725766169,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5323136\nallocs=16039\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 605970770,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5216624\nallocs=15610\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8732340,
            "unit": "ns",
            "extra": "gctime=0\nmemory=212560\nallocs=2051\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1136508660,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8979632\nallocs=55324\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1447091406.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42098944\nallocs=24956\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5495060,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73488\nallocs=1709\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 871912,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84176\nallocs=2043\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 5265492,
            "unit": "ns",
            "extra": "gctime=0\nmemory=244912\nallocs=5821\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3278086566,
            "unit": "ns",
            "extra": "gctime=0\nmemory=323519232\nallocs=4179231\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 155509232,
            "unit": "ns",
            "extra": "gctime=15741263\nmemory=169989864\nallocs=5995736\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 6783523.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/set_elements!",
            "value": 2453175,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_pair",
            "value": 65393333,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/iterator_interface",
            "value": 114608351.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=83202144\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_source",
            "value": 64928402,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/set_elements!",
            "value": 6199705,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_pair",
            "value": 265004037,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/iterator_interface",
            "value": 340453950,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_source",
            "value": 255345540,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2783957007.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71470480\nallocs=2241932\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 26849428,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131280\nallocs=2528\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "651849c3ff62458efb2b4fa1cbd243de64f7a6fd",
          "message": "Skip VectorOfVectors allocation test",
          "timestamp": "2026-01-07T11:03:22+01:00",
          "tree_id": "1a1afdf3e0f5f629a004f3d029efbb965242d131",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/651849c3ff62458efb2b4fa1cbd243de64f7a6fd"
        },
        "date": 1767780863861,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 13073823,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367264\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12935668,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 746215608,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5329184\nallocs=16141\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 627781926.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5216624\nallocs=15610\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8838561.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=212560\nallocs=2051\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1172161043,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9034064\nallocs=56242\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1454719942,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42153376\nallocs=25874\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5512317,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73488\nallocs=1709\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 875754,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84176\nallocs=2043\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 5296251,
            "unit": "ns",
            "extra": "gctime=0\nmemory=244912\nallocs=5821\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3250750351.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=330873680\nallocs=4199232\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 158501054,
            "unit": "ns",
            "extra": "gctime=15421122\nmemory=169989864\nallocs=5995736\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 6931341,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/set_elements!",
            "value": 2470657,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_pair",
            "value": 92305875.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/iterator_interface",
            "value": 147477627,
            "unit": "ns",
            "extra": "gctime=0\nmemory=83202144\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_source",
            "value": 94720636,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/set_elements!",
            "value": 8653232.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_pair",
            "value": 312419819.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/iterator_interface",
            "value": 401414437,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_source",
            "value": 313223418.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2868692866,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71476528\nallocs=2242034\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 28238781,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131280\nallocs=2528\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "ca9719936aa7375020a4006a2c4be6fab2145f9a",
          "message": "Tests: use Base.allocated instead of @allocated",
          "timestamp": "2026-01-07T12:35:19+01:00",
          "tree_id": "82a6668493db5db2a803cd99f1fb593a64ecb192",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/ca9719936aa7375020a4006a2c4be6fab2145f9a"
        },
        "date": 1767786163853,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 13146799,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367264\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12631654,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 709596406,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5317088\nallocs=15937\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 598708927,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5210576\nallocs=15508\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8638016,
            "unit": "ns",
            "extra": "gctime=0\nmemory=212560\nallocs=2051\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1155536855,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8979632\nallocs=55324\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1442948302.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42098944\nallocs=24956\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5486040,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73488\nallocs=1709\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 870123,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84176\nallocs=2043\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 5280325.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=244912\nallocs=5821\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3296847081,
            "unit": "ns",
            "extra": "gctime=0\nmemory=323521008\nallocs=4179250\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 156586040.5,
            "unit": "ns",
            "extra": "gctime=15996938.5\nmemory=169989864\nallocs=5995736\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 7104998,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/set_elements!",
            "value": 2415349,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_pair",
            "value": 68210074,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/iterator_interface",
            "value": 138214506,
            "unit": "ns",
            "extra": "gctime=0\nmemory=83202144\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_source",
            "value": 83049670,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/set_elements!",
            "value": 7497796,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_pair",
            "value": 252882927.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/iterator_interface",
            "value": 329060940,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_source",
            "value": 246425332,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2752771289,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71446288\nallocs=2241524\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 26972363,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131280\nallocs=2528\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "88a9465dbdd4b25e4022e096a2bf6f352b31fb7f",
          "message": "Allocation tests: try using closure",
          "timestamp": "2026-01-07T13:33:32+01:00",
          "tree_id": "5d94c883b9316089aa82da2a0c7610e79e939554",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/88a9465dbdd4b25e4022e096a2bf6f352b31fb7f"
        },
        "date": 1767789657676,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 12823361,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367264\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12690211.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 731333902,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5317088\nallocs=15937\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 610150643,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5210576\nallocs=15508\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8638482,
            "unit": "ns",
            "extra": "gctime=0\nmemory=212560\nallocs=2051\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1192990899,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8979632\nallocs=55324\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1463174442.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42098944\nallocs=24956\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5507597,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73488\nallocs=1709\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 870393,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84176\nallocs=2043\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 5279321,
            "unit": "ns",
            "extra": "gctime=0\nmemory=244912\nallocs=5821\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3191202992,
            "unit": "ns",
            "extra": "gctime=0\nmemory=323519744\nallocs=4179239\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 161014628,
            "unit": "ns",
            "extra": "gctime=16655851\nmemory=169989864\nallocs=5995736\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 6975711,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/set_elements!",
            "value": 2472200.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_pair",
            "value": 65318599,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/iterator_interface",
            "value": 101969436,
            "unit": "ns",
            "extra": "gctime=0\nmemory=83202144\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_source",
            "value": 70199793,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/set_elements!",
            "value": 6540716,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_pair",
            "value": 271266778,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/iterator_interface",
            "value": 349738774,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_source",
            "value": 268053095,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2806122507,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71446144\nallocs=2241523\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 27227292,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131280\nallocs=2528\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "241ac6958203eebcf994309bde015e840176c315",
          "message": "Use @ballocated in leapfrogging tests",
          "timestamp": "2026-01-07T15:39:58+01:00",
          "tree_id": "4937eff092806bd7477d1ba2aef07b3f30f3648a",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/241ac6958203eebcf994309bde015e840176c315"
        },
        "date": 1767797251358,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 13050482.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367264\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12963833,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 755030713,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5329184\nallocs=16141\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 634872479.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5216624\nallocs=15610\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8868521.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=212560\nallocs=2051\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1191870840,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9034064\nallocs=56242\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1489039663,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42153376\nallocs=25874\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5601730.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73488\nallocs=1709\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 875410,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84176\nallocs=2043\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 5439729,
            "unit": "ns",
            "extra": "gctime=0\nmemory=244912\nallocs=5821\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3258798663,
            "unit": "ns",
            "extra": "gctime=0\nmemory=323521520\nallocs=4179258\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 165599055,
            "unit": "ns",
            "extra": "gctime=16763858\nmemory=169989864\nallocs=5995736\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 6907374,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/set_elements!",
            "value": 2427026,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_pair",
            "value": 103869484,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/iterator_interface",
            "value": 137219888,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_source",
            "value": 122202172,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/set_elements!",
            "value": 9542697.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_pair",
            "value": 320634383,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/iterator_interface",
            "value": 405397443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_source",
            "value": 322098935,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2871427310.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71570368\nallocs=2244197\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 28463815,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131280\nallocs=2528\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "47347b0eb6f53b900deadeb371f721eebd7355ea",
          "message": "Fix tests?",
          "timestamp": "2026-01-07T17:28:04+01:00",
          "tree_id": "eb0faed3c18842325e79462ed0a603780a52043e",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/47347b0eb6f53b900deadeb371f721eebd7355ea"
        },
        "date": 1767803717432,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 12911961.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367264\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12617837,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 721436186,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5317088\nallocs=15937\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 608691794,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5210576\nallocs=15508\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8660571,
            "unit": "ns",
            "extra": "gctime=0\nmemory=212560\nallocs=2051\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1165056800,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8979632\nallocs=55324\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1435896966,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42098944\nallocs=24956\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5484494,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73488\nallocs=1709\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 872444.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84176\nallocs=2043\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 5257642.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=244912\nallocs=5821\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3243011900,
            "unit": "ns",
            "extra": "gctime=0\nmemory=323521008\nallocs=4179250\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 152092207,
            "unit": "ns",
            "extra": "gctime=15584859\nmemory=169989864\nallocs=5995736\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 7081667,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/set_elements!",
            "value": 2414350,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_pair",
            "value": 61344993,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/iterator_interface",
            "value": 97279300.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=83202144\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_source",
            "value": 61989636,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/set_elements!",
            "value": 5963434,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_pair",
            "value": 254686008.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/iterator_interface",
            "value": 328428022.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_source",
            "value": 253290966,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2777599361.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71446288\nallocs=2241524\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 26647710,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131280\nallocs=2528\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "ad70238e87125f0f2cc99a4cfb766a15c1b279ac",
          "message": "Add scripts/merge_vtkhdf_series.jl",
          "timestamp": "2026-01-12T11:57:52+01:00",
          "tree_id": "ad3d6e0fc021ad1609aa81257e0203b09b875d3c",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/ad70238e87125f0f2cc99a4cfb766a15c1b279ac"
        },
        "date": 1768216318057,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 12483365,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367264\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12193090,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 707366936.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5323136\nallocs=16039\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 599925905,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5216624\nallocs=15610\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8427528.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=212560\nallocs=2051\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1094218835,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8979632\nallocs=55324\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1458396463,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42098944\nallocs=24956\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5387650,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73488\nallocs=1709\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 871243,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84176\nallocs=2043\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 5238062,
            "unit": "ns",
            "extra": "gctime=0\nmemory=244912\nallocs=5821\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3296144492,
            "unit": "ns",
            "extra": "gctime=0\nmemory=323521520\nallocs=4179258\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 145171433.5,
            "unit": "ns",
            "extra": "gctime=15019083.5\nmemory=169989864\nallocs=5995736\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 6685008,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/set_elements!",
            "value": 2416509,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_pair",
            "value": 57688895,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/iterator_interface",
            "value": 93463517,
            "unit": "ns",
            "extra": "gctime=0\nmemory=83202144\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_source",
            "value": 59212438,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/set_elements!",
            "value": 5805218,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_pair",
            "value": 245937381,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/iterator_interface",
            "value": 316613210.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_source",
            "value": 243748752,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2712042516.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71470480\nallocs=2241932\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 27403124.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131280\nallocs=2528\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "d812f69ff0207515bef17ddbcbbcff280417bb9c",
          "message": "Fix folder name",
          "timestamp": "2026-01-12T12:55:28+01:00",
          "tree_id": "e358724a834d3765592f8bbb2280f5190685c7a9",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/d812f69ff0207515bef17ddbcbbcff280417bb9c"
        },
        "date": 1768219384956,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 13041368,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367264\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12923559,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 743528431,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5335232\nallocs=16243\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 651083014,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5222672\nallocs=15712\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8811641,
            "unit": "ns",
            "extra": "gctime=0\nmemory=212560\nallocs=2051\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1172690544,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9034064\nallocs=56242\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1631213339,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42153376\nallocs=25874\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5512570.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73488\nallocs=1709\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 871257,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84176\nallocs=2043\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 5283590,
            "unit": "ns",
            "extra": "gctime=0\nmemory=244912\nallocs=5821\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3194920074,
            "unit": "ns",
            "extra": "gctime=0\nmemory=323522032\nallocs=4179266\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 161705430,
            "unit": "ns",
            "extra": "gctime=16615780\nmemory=169990008\nallocs=5995739\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 6885832,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/set_elements!",
            "value": 2417617,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_pair",
            "value": 65711365,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/iterator_interface",
            "value": 108240672,
            "unit": "ns",
            "extra": "gctime=0\nmemory=83202144\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_source",
            "value": 66552750.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/set_elements!",
            "value": 6580465,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_pair",
            "value": 275750164,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/iterator_interface",
            "value": 359582112.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_source",
            "value": 274249666,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2822936042,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71500720\nallocs=2242442\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 27832044.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131280\nallocs=2528\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "37919337ed63378b050deebd08e8d56feb6bfa67",
          "message": "Remove Core.Box in write_vtkhdf\n\nDetected using script from https://github.com/JuliaLang/julia/pull/60478\n\nThis was the only Core.Box detected in VortexPasta.",
          "timestamp": "2026-01-12T13:10:11+01:00",
          "tree_id": "1911dfa7b892c12cf9adf41766427640204f22eb",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/37919337ed63378b050deebd08e8d56feb6bfa67"
        },
        "date": 1768220314899,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 13296162.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367264\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 13021143,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 746809056,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5329184\nallocs=16141\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 627543807.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5216624\nallocs=15610\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8741849.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=212560\nallocs=2051\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1238821534,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9034064\nallocs=56242\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1508637088,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42153376\nallocs=25874\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5506415,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73488\nallocs=1709\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 872313,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84176\nallocs=2043\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 5282774.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=244912\nallocs=5821\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3251814110,
            "unit": "ns",
            "extra": "gctime=0\nmemory=323521264\nallocs=4179254\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 166414153,
            "unit": "ns",
            "extra": "gctime=16223525.5\nmemory=169989864\nallocs=5995736\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 7054021,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/set_elements!",
            "value": 2433317,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_pair",
            "value": 68961077.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/iterator_interface",
            "value": 130176025,
            "unit": "ns",
            "extra": "gctime=0\nmemory=83202144\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_source",
            "value": 69913095,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/set_elements!",
            "value": 7893717,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_pair",
            "value": 296719735,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/iterator_interface",
            "value": 399568842,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_source",
            "value": 285995474,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2937775185,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71570968\nallocs=2242199\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 28422126,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131584\nallocs=2534\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "60b61339876edb5961c7aa6ec0768c8b5d4bc628",
          "message": "Tests: avoid Optim.jl v2 for now",
          "timestamp": "2026-01-12T14:25:30+01:00",
          "tree_id": "d3f425ac5d7e4f7d25aa0210cd4df70279da4741",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/60b61339876edb5961c7aa6ec0768c8b5d4bc628"
        },
        "date": 1768224790834,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 13150477.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367264\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12999119.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 771679549,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5320480\nallocs=16137\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 646497065.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5207920\nallocs=15606\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8976085,
            "unit": "ns",
            "extra": "gctime=0\nmemory=212560\nallocs=2051\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1255285859,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9034064\nallocs=56242\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1510207267,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42153376\nallocs=25874\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5553788.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73488\nallocs=1709\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 874428,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84176\nallocs=2043\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 5327868.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=244912\nallocs=5821\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3202978717.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=323521776\nallocs=4179262\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 169272291,
            "unit": "ns",
            "extra": "gctime=17432935\nmemory=169989864\nallocs=5995736\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 7240055.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/set_elements!",
            "value": 2480975,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_pair",
            "value": 124834075,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/iterator_interface",
            "value": 90105134,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_source",
            "value": 74396079.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/set_elements!",
            "value": 9327178.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_pair",
            "value": 335409960.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/iterator_interface",
            "value": 413708861,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_source",
            "value": 315745798.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2849571116.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71540592\nallocs=2244490\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 28973588,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131280\nallocs=2528\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "845af707082a9155e95516b89632eb7f780037ff",
          "message": "Don't use @ballocated in ring_energy test",
          "timestamp": "2026-01-12T17:09:51+01:00",
          "tree_id": "8ea14bb23a95c0ded7a74753f6408a5457053c4d",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/845af707082a9155e95516b89632eb7f780037ff"
        },
        "date": 1768234627702,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 12831939,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367264\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12675642,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 718482030,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5317088\nallocs=15937\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 602002651,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5210576\nallocs=15508\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8659146,
            "unit": "ns",
            "extra": "gctime=0\nmemory=212560\nallocs=2051\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1113869598,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8979632\nallocs=55324\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1506849962,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42098944\nallocs=24956\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5512014,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73488\nallocs=1709\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 873440.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84176\nallocs=2043\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 5407514,
            "unit": "ns",
            "extra": "gctime=0\nmemory=244912\nallocs=5821\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3175615677.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=323518480\nallocs=4179216\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 154781009,
            "unit": "ns",
            "extra": "gctime=16114209\nmemory=169989864\nallocs=5995736\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 7114680,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/set_elements!",
            "value": 2425467,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_pair",
            "value": 62578148,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/iterator_interface",
            "value": 97435092.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=83202144\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_source",
            "value": 62731873,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/set_elements!",
            "value": 6052329,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_pair",
            "value": 252613267,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/iterator_interface",
            "value": 326791253.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_source",
            "value": 252200293,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2790399995.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71446288\nallocs=2241524\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 28724753,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131280\nallocs=2528\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "9bb0fbe10aaf169d581dc66db5859fd503094ef2",
          "message": "Fix forcing docs",
          "timestamp": "2026-01-15T12:17:40+01:00",
          "tree_id": "01de877189ac1a62295a997bc7a8d12e4179fd6f",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/9bb0fbe10aaf169d581dc66db5859fd503094ef2"
        },
        "date": 1768476535158,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 13115581.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367264\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 13002367,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 747317488,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5317088\nallocs=15937\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 629888875,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5210576\nallocs=15508\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8949397,
            "unit": "ns",
            "extra": "gctime=0\nmemory=212560\nallocs=2051\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1164297998,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8979632\nallocs=55324\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1525622808.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42098944\nallocs=24956\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5535785,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73488\nallocs=1709\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 877939.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84176\nallocs=2043\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 5267020,
            "unit": "ns",
            "extra": "gctime=0\nmemory=244912\nallocs=5821\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3190540805.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=330874416\nallocs=4199242\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 159008638,
            "unit": "ns",
            "extra": "gctime=16282789\nmemory=169990008\nallocs=5995739\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 6976760,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/set_elements!",
            "value": 2439043,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_pair",
            "value": 70242663.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/iterator_interface",
            "value": 130466176,
            "unit": "ns",
            "extra": "gctime=0\nmemory=83202144\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_source",
            "value": 73407864.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/set_elements!",
            "value": 8207731,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_pair",
            "value": 330734775,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/iterator_interface",
            "value": 413312203,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_source",
            "value": 318601761.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2916646446.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71446288\nallocs=2241524\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 28013357.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131280\nallocs=2528\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "7caa8108364ca5cfe04646116f2fc7ca7263a25b",
          "message": "Disable some allocation tests on Julia 1.12",
          "timestamp": "2026-01-15T13:46:05+01:00",
          "tree_id": "7fb3c833c9affff81fc052ee83e1457be30aa6fc",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/7caa8108364ca5cfe04646116f2fc7ca7263a25b"
        },
        "date": 1768481608324,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 12971821.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367264\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12706288.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 716238943,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5323136\nallocs=16039\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 602040517,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5216624\nallocs=15610\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8647141,
            "unit": "ns",
            "extra": "gctime=0\nmemory=212560\nallocs=2051\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1149849833,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8979632\nallocs=55324\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1518132074.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42098944\nallocs=24956\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5497193,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4016\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 870716,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84176\nallocs=2043\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 5273059.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=244912\nallocs=5821\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3253387068.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=323521104\nallocs=4179256\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 155189906,
            "unit": "ns",
            "extra": "gctime=16640554\nmemory=169989864\nallocs=5995736\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 6751779,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/set_elements!",
            "value": 2419541.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_pair",
            "value": 63614218.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/iterator_interface",
            "value": 107851983,
            "unit": "ns",
            "extra": "gctime=0\nmemory=83202144\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_source",
            "value": 64530364,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/set_elements!",
            "value": 6169967,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_pair",
            "value": 260400066.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/iterator_interface",
            "value": 335140961,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_source",
            "value": 258801504.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2772613860.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71564920\nallocs=2242097\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 29644959.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131584\nallocs=2534\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "52d0a1b61a5d74ece072d203f00770edea5665af",
          "message": "Docs: update load_checkpoint example",
          "timestamp": "2026-01-20T13:40:12+01:00",
          "tree_id": "d77df5579d9e46450a32e8d5e6135e3dfe48d932",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/52d0a1b61a5d74ece072d203f00770edea5665af"
        },
        "date": 1768913470182,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 13107396,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367264\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12912566,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 740111738,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5317088\nallocs=15937\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 627219509,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5210576\nallocs=15508\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8917317,
            "unit": "ns",
            "extra": "gctime=0\nmemory=212560\nallocs=2051\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1126938139,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8979632\nallocs=55324\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1583960570,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42098944\nallocs=24956\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5502236,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73488\nallocs=1709\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 875706,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84176\nallocs=2043\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 5252760,
            "unit": "ns",
            "extra": "gctime=0\nmemory=244912\nallocs=5821\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3229884052.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=323520752\nallocs=4179246\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 155203666,
            "unit": "ns",
            "extra": "gctime=15630715.5\nmemory=169989864\nallocs=5995736\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 6810703,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/set_elements!",
            "value": 2422591,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_pair",
            "value": 67678876,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/iterator_interface",
            "value": 117728163,
            "unit": "ns",
            "extra": "gctime=0\nmemory=83202144\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_source",
            "value": 71579669,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/set_elements!",
            "value": 6674483,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_pair",
            "value": 270565470,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/iterator_interface",
            "value": 353497218,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_source",
            "value": 267688691,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2837107167.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71540128\nallocs=2243687\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 28046690.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131280\nallocs=2528\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "33ac40fbf5a7752f29d3925eb034a5498d88786f",
          "message": "Some changes to spectrum implementation",
          "timestamp": "2026-01-21T15:35:36+01:00",
          "tree_id": "6e3fdc209e7b6accddd694998dd4b0aa61097841",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/33ac40fbf5a7752f29d3925eb034a5498d88786f"
        },
        "date": 1769006785008,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 12730630,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367264\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12638419,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 730978346,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5308384\nallocs=15933\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 616713648,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5201872\nallocs=15504\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8655414,
            "unit": "ns",
            "extra": "gctime=0\nmemory=212560\nallocs=2051\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1116535072,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8979632\nallocs=55324\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1402592871,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42098944\nallocs=24956\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5484414,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73488\nallocs=1709\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 874120.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84176\nallocs=2043\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 3963734.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5312\nallocs=62\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3212434676.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=323520912\nallocs=4179250\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 149901348,
            "unit": "ns",
            "extra": "gctime=15169042\nmemory=169989864\nallocs=5995736\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 6936781,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/set_elements!",
            "value": 2437085.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_pair",
            "value": 60448699,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/iterator_interface",
            "value": 97195043.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=83202144\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_source",
            "value": 60587877,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/set_elements!",
            "value": 5941938.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_pair",
            "value": 248584594,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/iterator_interface",
            "value": 319195885,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_source",
            "value": 249100039,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2812210665,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71605112\nallocs=2244153\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 27994948,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131584\nallocs=2534\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "0cfa0a45e21f9968778eb52ddba44324eb0c3d20",
          "message": "Fix spectrum on CPU",
          "timestamp": "2026-01-22T10:55:08+01:00",
          "tree_id": "db093b57999a8369948ea6d083811b9e2c48157a",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/0cfa0a45e21f9968778eb52ddba44324eb0c3d20"
        },
        "date": 1769076721662,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 12711091,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367264\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12640224,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 727133551,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5323136\nallocs=16039\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 614055156,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5216624\nallocs=15610\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8649067,
            "unit": "ns",
            "extra": "gctime=0\nmemory=212560\nallocs=2051\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1141887769,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8979632\nallocs=55324\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1517065374,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42098944\nallocs=24956\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5511383.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73488\nallocs=1709\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 871864,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84176\nallocs=2043\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 3949507,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5312\nallocs=62\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3209494134,
            "unit": "ns",
            "extra": "gctime=0\nmemory=323520144\nallocs=4179240\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 154337220,
            "unit": "ns",
            "extra": "gctime=16010912.5\nmemory=169989864\nallocs=5995736\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 6738529.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/set_elements!",
            "value": 2412012.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_pair",
            "value": 63074960,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/iterator_interface",
            "value": 99637686.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=83202144\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_source",
            "value": 63686350,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/set_elements!",
            "value": 6038339.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_pair",
            "value": 256090318.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/iterator_interface",
            "value": 331817451.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_source",
            "value": 256743104,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2815522556.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71470480\nallocs=2241932\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 26041788,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131280\nallocs=2528\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "juan-ignacio.polanco@cnrs.fr",
            "name": "Juan Ignacio Polanco",
            "username": "jipolanco"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "80b9ffa26b0c829c8f9b37c4ce8a3f57872edff1",
          "message": "Split CPU and GPU implementations of energy spectrum (#88)\n\n* Split CPU and GPU implementations of energy spectrum\n\nAlso, simplify *and* speed-up GPU implementation, using atomics instead\nof shared memory.\n\n* Cleanup\n\n* Tune GPU group size\n\n* Remove large allocation in CPU spectrum\n\n* Use @threads instead of @spawn\n\n* Remove PseudoGPU",
          "timestamp": "2026-01-22T15:51:50+01:00",
          "tree_id": "a4fb9acac87361f7b28623600c11f26447ef0d04",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/80b9ffa26b0c829c8f9b37c4ce8a3f57872edff1"
        },
        "date": 1769094392837,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 13107038,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367264\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12845929,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 753500676,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5317088\nallocs=15937\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 620543259,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5210576\nallocs=15508\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8776427.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=212560\nallocs=2051\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1151175943,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8979632\nallocs=55324\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1422000284,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42098944\nallocs=24956\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5502464,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4016\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 872262,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84176\nallocs=2043\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2273431,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3159079697,
            "unit": "ns",
            "extra": "gctime=0\nmemory=323521760\nallocs=4179261\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 161753655,
            "unit": "ns",
            "extra": "gctime=16227842\nmemory=169989864\nallocs=5995736\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 7173223,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/set_elements!",
            "value": 2433517,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_pair",
            "value": 76121027,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/iterator_interface",
            "value": 129096337,
            "unit": "ns",
            "extra": "gctime=0\nmemory=83202144\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_source",
            "value": 74585315,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/set_elements!",
            "value": 7775609.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_pair",
            "value": 313583949.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/iterator_interface",
            "value": 398331815,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_source",
            "value": 300998574,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2842935712,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71540584\nallocs=2241688\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 28742328,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131584\nallocs=2534\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "e74f432b47a3ac0d4bc09726ff049fc4f055931f",
          "message": "v0.32.13",
          "timestamp": "2026-01-22T16:37:12+01:00",
          "tree_id": "be1e0449d74c65b046188b0e0b6a1b324424f5e6",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/e74f432b47a3ac0d4bc09726ff049fc4f055931f"
        },
        "date": 1769096869815,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 12939701.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367264\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12703698.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 730323807,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5317088\nallocs=15937\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 617590161,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5210576\nallocs=15508\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8686653,
            "unit": "ns",
            "extra": "gctime=0\nmemory=212560\nallocs=2051\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1151510322,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8979632\nallocs=55324\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1528230522,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42098944\nallocs=24956\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5502666,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73488\nallocs=1709\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 872266.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84176\nallocs=2043\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2254306,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3261508527.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=323521664\nallocs=4179263\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 152729325,
            "unit": "ns",
            "extra": "gctime=15028486.5\nmemory=169989864\nallocs=5995736\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 6842044,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/set_elements!",
            "value": 2458791,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_pair",
            "value": 65684880.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/iterator_interface",
            "value": 78476046,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_source",
            "value": 66391202,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/set_elements!",
            "value": 6188774,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_pair",
            "value": 260875793.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/iterator_interface",
            "value": 335406286,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_source",
            "value": 258511024,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2828025694.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71446288\nallocs=2241524\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 28230143,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131280\nallocs=2528\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "a04823bbeef5d1c121a3e3110c94a6528a74ac93",
          "message": "Allow optional function in load_checkpoint",
          "timestamp": "2026-01-23T14:40:35+01:00",
          "tree_id": "63b6256751ff06b29426d7f3166a1a4fe0ae411f",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/a04823bbeef5d1c121a3e3110c94a6528a74ac93"
        },
        "date": 1769176118693,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 12957388,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367264\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12806417.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 739970753,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5329184\nallocs=16141\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 613996657,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5216624\nallocs=15610\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8664438.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=212560\nallocs=2051\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1163150276,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9034064\nallocs=56242\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1612899526.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42153376\nallocs=25874\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5498537.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73488\nallocs=1709\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 871199,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84176\nallocs=2043\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2286554,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3194523755,
            "unit": "ns",
            "extra": "gctime=0\nmemory=323520656\nallocs=4179248\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 167276803.5,
            "unit": "ns",
            "extra": "gctime=16219911\nmemory=169989864\nallocs=5995736\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 7206987,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/set_elements!",
            "value": 2415795,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_pair",
            "value": 80302363,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/iterator_interface",
            "value": 126815631.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=83202144\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_source",
            "value": 66846690,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/set_elements!",
            "value": 7173043,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_pair",
            "value": 291135810,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/iterator_interface",
            "value": 379967574,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_source",
            "value": 294553913,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2817540137,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71476528\nallocs=2242034\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 28890805,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131280\nallocs=2528\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "b5bb7db25a9e0b372e2701dba3fb08476b19b434",
          "message": "v0.32.14",
          "timestamp": "2026-01-23T16:41:26+01:00",
          "tree_id": "cb747abd9db82914e4cec7766ea113fcbccf382e",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/b5bb7db25a9e0b372e2701dba3fb08476b19b434"
        },
        "date": 1769183316828,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 12880532,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367264\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12794390,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 726024502,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5308384\nallocs=15933\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 614651162,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5201872\nallocs=15504\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8719508,
            "unit": "ns",
            "extra": "gctime=0\nmemory=212560\nallocs=2051\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1156959126,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8979632\nallocs=55324\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1617560339.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42098944\nallocs=24956\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5500935,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73488\nallocs=1709\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 872482,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84176\nallocs=2043\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2263819,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3220346486.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=323521264\nallocs=4179254\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 156580224,
            "unit": "ns",
            "extra": "gctime=16389815.5\nmemory=169989864\nallocs=5995736\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 6799280,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/set_elements!",
            "value": 2436901.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_pair",
            "value": 63760038.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/iterator_interface",
            "value": 112341485.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=83202144\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_source",
            "value": 65825647,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/set_elements!",
            "value": 6700829,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_pair",
            "value": 267924081,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/iterator_interface",
            "value": 343420744,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_source",
            "value": 270858151,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2840648008,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71411472\nallocs=2241508\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 28909500,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131280\nallocs=2528\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "b9857a566989349ecb28f7a016370f33977da32e",
          "message": "Refactor FourierBandForcingBS",
          "timestamp": "2026-01-26T11:09:09+01:00",
          "tree_id": "4027ff0d49c41d0cd4d7a1fc8e11877323ab1778",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/b9857a566989349ecb28f7a016370f33977da32e"
        },
        "date": 1769422644001,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 13005318,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367264\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12984663,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 763500312,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5329184\nallocs=16141\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 646577795,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5216624\nallocs=15610\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8925806,
            "unit": "ns",
            "extra": "gctime=0\nmemory=212560\nallocs=2051\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1229937779,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9034064\nallocs=56242\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1676861147,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42153376\nallocs=25874\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5576468.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73488\nallocs=1709\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 877779,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84176\nallocs=2043\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2268878,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3170940857,
            "unit": "ns",
            "extra": "gctime=0\nmemory=323520576\nallocs=4179251\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 161608072,
            "unit": "ns",
            "extra": "gctime=15938624\nmemory=169989864\nallocs=5995736\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 7296494,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/set_elements!",
            "value": 2416737,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_pair",
            "value": 101398518,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/iterator_interface",
            "value": 125709411,
            "unit": "ns",
            "extra": "gctime=0\nmemory=83202144\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 1/foreach_source",
            "value": 74166092,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/set_elements!",
            "value": 6815052,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_pair",
            "value": 290619970,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/iterator_interface",
            "value": 379386198.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/nsubdiv = 2/foreach_source",
            "value": 265772647,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2839962816.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71476528\nallocs=2242034\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 28901145,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131280\nallocs=2528\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "947628715edf27b9d5230c9f95c4e4fcc9ba9552",
          "message": "Update CellLists benchmarks",
          "timestamp": "2026-01-26T12:09:31+01:00",
          "tree_id": "f131e4b168cb42e9294e507471fbb12de05ad822",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/947628715edf27b9d5230c9f95c4e4fcc9ba9552"
        },
        "date": 1769426254029,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 13315874,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367264\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 13113763,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 754275541,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5308384\nallocs=15933\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 618856964,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5201872\nallocs=15504\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8711033,
            "unit": "ns",
            "extra": "gctime=0\nmemory=212560\nallocs=2051\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1158210466,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8979632\nallocs=55324\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1479748255,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42098944\nallocs=24956\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5568048,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73488\nallocs=1709\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 876558.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84176\nallocs=2043\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2274582,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3258505360.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=323522368\nallocs=4179267\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 163222293,
            "unit": "ns",
            "extra": "gctime=16332949\nmemory=169989864\nallocs=5995736\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 6925294,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 2014897,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair",
            "value": 74869626,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 128399231,
            "unit": "ns",
            "extra": "gctime=0\nmemory=83202144\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 80966934.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD)",
            "value": 634515993.5,
            "unit": "ns",
            "extra": "gctime=238793804.5\nmemory=1568657872\nallocs=14655065\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 6776392.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair",
            "value": 283731778,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 378854628.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 292578394,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD)",
            "value": 298049188,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2829170714,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71406432\nallocs=2241199\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 28191704,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131280\nallocs=2528\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "fbb775472c75d295e9f41e3076ba19f1c882fd5b",
          "message": "Add CellLists benchmarks using OpenCL",
          "timestamp": "2026-01-26T12:24:35+01:00",
          "tree_id": "f38d8d786e785e2822c643105d8087ba64866f81",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/fbb775472c75d295e9f41e3076ba19f1c882fd5b"
        },
        "date": 1769427286266,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 12933144.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367264\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12721331,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 720000524,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5317088\nallocs=15937\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 608565531,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5210576\nallocs=15508\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8663996.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=212560\nallocs=2051\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1123052584,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8979632\nallocs=55324\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1585377622.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42098944\nallocs=24956\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5510083,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73488\nallocs=1709\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 832538,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4016\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2509345.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3169769188,
            "unit": "ns",
            "extra": "gctime=0\nmemory=323520656\nallocs=4179248\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 153574292.5,
            "unit": "ns",
            "extra": "gctime=17429486.5\nmemory=169990008\nallocs=5995739\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 6905497,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 2020551,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair",
            "value": 62578971.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 101651260,
            "unit": "ns",
            "extra": "gctime=0\nmemory=83202144\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 62742358,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD)",
            "value": 86892241,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 5795360.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3648\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair",
            "value": 257706740,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 331849051,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 253874382,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD)",
            "value": 255995633,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/set_elements!",
            "value": 2359089,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45392\nallocs=731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/foreach_pair",
            "value": 64069359,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/set_elements!",
            "value": 7223648,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45392\nallocs=731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/foreach_pair",
            "value": 207464338,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2775690905.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=71446288\nallocs=2241524\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 27219793,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131280\nallocs=2528\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "61061062a70ee8ba0c8ed0d44f72c31781608517",
          "message": "v0.32.15",
          "timestamp": "2026-01-27T10:27:39+01:00",
          "tree_id": "fce41f0f21d0209b4b83d12846627a4ab44ce2d0",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/61061062a70ee8ba0c8ed0d44f72c31781608517"
        },
        "date": 1769507207726,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 13261625,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367520\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12610358.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 768221137,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5339648\nallocs=16287\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 647967308.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5226896\nallocs=15756\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8816750.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=213056\nallocs=2049\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1267176909.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9054936\nallocs=56368\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1698773122,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42161936\nallocs=25964\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5498258.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73696\nallocs=1714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 873278.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84416\nallocs=2049\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2271117.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3164746387,
            "unit": "ns",
            "extra": "gctime=0\nmemory=334882480\nallocs=4179262\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 303154265,
            "unit": "ns",
            "extra": "gctime=26244953\nmemory=209976264\nallocs=7493721\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 7010756,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 2008942,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (sorted)",
            "value": 61154775,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/unsorted)",
            "value": 97648119,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 84381400,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 68760571,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/sorted)",
            "value": 79870707,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 72293076,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 6221798,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (sorted)",
            "value": 108096710,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/unsorted)",
            "value": 280942153.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 382675691,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 289933762,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/sorted)",
            "value": 124193869,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 309533352,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/set_elements!",
            "value": 2354580,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45392\nallocs=731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/foreach_pair (sorted)",
            "value": 65776083,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41696\nallocs=641\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 68873128,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/set_elements!",
            "value": 7517684,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45392\nallocs=731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/foreach_pair (sorted)",
            "value": 191019283.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41696\nallocs=641\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 227735508,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 3095698943.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68944288\nallocs=2242325\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 27322811,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131888\nallocs=2540\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "juan-ignacio.polanco@cnrs.fr",
            "name": "Juan Ignacio Polanco",
            "username": "jipolanco"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "5966adab751d214637d1b517cbfa6035f888b013",
          "message": "Add SIMD implementation of exp (#90)\n\n* Add SIMD implementation of exp\n\n* Update CHANGELOG\n\n* Fix exp Float32 + add tests\n\n* Fix non-periodic + SIMD\n\n* Port vectorised erfc from xsimd",
          "timestamp": "2026-01-28T14:59:10+01:00",
          "tree_id": "4d1dfc151bb317a88cb2bd07f9c348cccfabecdc",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/5966adab751d214637d1b517cbfa6035f888b013"
        },
        "date": 1769609860605,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 13357724,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367520\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12825047,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 710327524,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5321504\nallocs=15981\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 602682911,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5214864\nallocs=15552\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8665257,
            "unit": "ns",
            "extra": "gctime=0\nmemory=213184\nallocs=2057\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1155957534,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9002808\nallocs=55594\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1533188950,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42107504\nallocs=25046\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5495057.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73696\nallocs=1714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 875233.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84416\nallocs=2049\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2247791.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3214594055,
            "unit": "ns",
            "extra": "gctime=0\nmemory=334879984\nallocs=4179230\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 262405550.5,
            "unit": "ns",
            "extra": "gctime=24737796\nmemory=209976264\nallocs=7493721\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 6876030,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 2007740.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (sorted)",
            "value": 56380241,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/unsorted)",
            "value": 89670902,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 79451255,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 65996353,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/sorted)",
            "value": 77180983,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 65171574,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 5879095.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (sorted)",
            "value": 100959734,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/unsorted)",
            "value": 269961043,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 351803141,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 262689659.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/sorted)",
            "value": 106179489,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 273047921,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/set_elements!",
            "value": 2315786.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45392\nallocs=731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/foreach_pair (sorted)",
            "value": 62798437.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41696\nallocs=641\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 63859621,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/set_elements!",
            "value": 7800227,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45392\nallocs=731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/foreach_pair (sorted)",
            "value": 191158570,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41696\nallocs=641\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 228710607.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2804078705,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68890048\nallocs=2241407\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 29028556,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131888\nallocs=2540\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "1db7d91202378060800670360b0534770a959863",
          "message": "v0.32.16",
          "timestamp": "2026-01-28T15:00:58+01:00",
          "tree_id": "b89e6441efe5ce1c85eef9cd3614b4d4db99731b",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/1db7d91202378060800670360b0534770a959863"
        },
        "date": 1769610009236,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 13674746,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367520\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 13115180,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 729950700,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5339648\nallocs=16287\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 609694789,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5226960\nallocs=15756\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8770449,
            "unit": "ns",
            "extra": "gctime=0\nmemory=213184\nallocs=2057\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1164758230,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9057240\nallocs=56512\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1521341488,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42161936\nallocs=25964\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5516437,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73696\nallocs=1714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 877181,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84416\nallocs=2049\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2269120,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3187968194,
            "unit": "ns",
            "extra": "gctime=0\nmemory=334882992\nallocs=4179270\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 266667345.5,
            "unit": "ns",
            "extra": "gctime=25520202.5\nmemory=209976264\nallocs=7493721\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 7275661,
            "unit": "ns",
            "extra": "gctime=0\nmemory=851184\nallocs=19584\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 2067399,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (sorted)",
            "value": 59899528.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/unsorted)",
            "value": 108622451,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 119454268,
            "unit": "ns",
            "extra": "gctime=0\nmemory=89602208\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 69752488,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/sorted)",
            "value": 79599624,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 75642656,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 6930912,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (sorted)",
            "value": 102309553,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/unsorted)",
            "value": 285657972.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 385519538,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 291637760.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/sorted)",
            "value": 121709202.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 315384727,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/set_elements!",
            "value": 2384931,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45392\nallocs=731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/foreach_pair (sorted)",
            "value": 68469998,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41696\nallocs=641\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 69596512,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/set_elements!",
            "value": 8216979,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45392\nallocs=731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/foreach_pair (sorted)",
            "value": 190443236,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41696\nallocs=641\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 223794590,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2790611527.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68909632\nallocs=2242761\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 28698994,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131888\nallocs=2540\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "juan-ignacio.polanco@cnrs.fr",
            "name": "Juan Ignacio Polanco",
            "username": "jipolanco"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "ee294bca97a5b139d7e60a084ea631946f6a9a71",
          "message": "Activate GPU device before resizing arrays (#91)",
          "timestamp": "2026-01-29T13:15:33+01:00",
          "tree_id": "7e4dbef663985a021f61933a25e6cfaba3e2f11c",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/ee294bca97a5b139d7e60a084ea631946f6a9a71"
        },
        "date": 1769690058645,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 13348838,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367520\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12801852,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 714494968,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5312608\nallocs=15973\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 595533019,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5206160\nallocs=15548\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8667783.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=213184\nallocs=2057\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1136781029,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9002808\nallocs=55594\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1576933761.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42107504\nallocs=25046\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5484169,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73696\nallocs=1714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 874500,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84416\nallocs=2049\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2247569,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3216710320.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=334884096\nallocs=4179283\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 265787158.5,
            "unit": "ns",
            "extra": "gctime=25778291\nmemory=209976264\nallocs=7493721\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 6910270,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 1972409,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (sorted)",
            "value": 57448488,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/unsorted)",
            "value": 102826719.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 125602948,
            "unit": "ns",
            "extra": "gctime=0\nmemory=89602208\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 67296696,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/sorted)",
            "value": 77467012,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 73632698,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 8108608,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (sorted)",
            "value": 110591759,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/unsorted)",
            "value": 290179958,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 379354074,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 291797782,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/sorted)",
            "value": 114171669,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 313902915,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/set_elements!",
            "value": 2379634,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45392\nallocs=731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/foreach_pair (sorted)",
            "value": 63588832,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41696\nallocs=641\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 65082218,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/set_elements!",
            "value": 7777754,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45392\nallocs=731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/foreach_pair (sorted)",
            "value": 189466825,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41696\nallocs=641\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 236856539,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2819456264.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68855040\nallocs=2241387\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 28659044,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131888\nallocs=2540\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "8733b9fded60ea97c4cc95048f5360030353fc86",
          "message": "v0.32.17",
          "timestamp": "2026-01-29T13:17:43+01:00",
          "tree_id": "57459114339a37ccca3573f0532f8f95a412b846",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/8733b9fded60ea97c4cc95048f5360030353fc86"
        },
        "date": 1769690170914,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 13463228,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367520\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12644332,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 706674929.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5333408\nallocs=16181\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 585966977,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5220912\nallocs=15654\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8675041,
            "unit": "ns",
            "extra": "gctime=0\nmemory=213184\nallocs=2057\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1143774874,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9057240\nallocs=56512\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1506077749.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42161936\nallocs=25964\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5485299,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73696\nallocs=1714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 873720.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84416\nallocs=2049\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2242870.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3249997117.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=334883328\nallocs=4179271\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 263181921,
            "unit": "ns",
            "extra": "gctime=25686487.5\nmemory=209976264\nallocs=7493721\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 6917151,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 2056564.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (sorted)",
            "value": 54119706.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/unsorted)",
            "value": 83613746.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 97818862,
            "unit": "ns",
            "extra": "gctime=0\nmemory=89602208\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 61855836,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/sorted)",
            "value": 73235284,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 61122097,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 5323788,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (sorted)",
            "value": 94404679,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/unsorted)",
            "value": 249438078,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 318174693.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 242914899,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/sorted)",
            "value": 105563792,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 245120118,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/set_elements!",
            "value": 2368466,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45392\nallocs=731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/foreach_pair (sorted)",
            "value": 62306336,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41696\nallocs=641\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 66139267,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/set_elements!",
            "value": 6987768.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45392\nallocs=731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/foreach_pair (sorted)",
            "value": 182705757.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41696\nallocs=641\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 202082031,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2759039166,
            "unit": "ns",
            "extra": "gctime=0\nmemory=69019328\nallocs=2244399\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 28692202,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131888\nallocs=2540\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "juan-ignacio.polanco@cnrs.fr",
            "name": "Juan Ignacio Polanco",
            "username": "jipolanco"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "da5f0087adb1afd3ac7422829fbae142b9d39f3d",
          "message": "Show movie/gif of turbulence simulation in docs and README (#92)\n\n* README: add embedded movie\n\n* Use gif instead\n\n* Docs: embed QT movie from vimeo\n\n* README: make the gif clickable",
          "timestamp": "2026-01-29T17:18:39+01:00",
          "tree_id": "74783a40c24a8f08d75a701ab05c153001a8dd31",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/da5f0087adb1afd3ac7422829fbae142b9d39f3d"
        },
        "date": 1769704368478,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 13533983.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367520\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12655123,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 708625688,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5312608\nallocs=15973\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 591211781,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5206160\nallocs=15548\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8675964,
            "unit": "ns",
            "extra": "gctime=0\nmemory=213184\nallocs=2057\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1139078531,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9002808\nallocs=55594\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1418478994,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42107504\nallocs=25046\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5493727,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73696\nallocs=1714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 875211.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84416\nallocs=2049\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2246445,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3217790657.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=334882064\nallocs=4179260\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 262444457.5,
            "unit": "ns",
            "extra": "gctime=26505683\nmemory=209976264\nallocs=7493721\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 6790746,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 1986545,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (sorted)",
            "value": 54263113,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/unsorted)",
            "value": 86294515.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 99385306,
            "unit": "ns",
            "extra": "gctime=0\nmemory=89602208\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 62322865,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/sorted)",
            "value": 74042838.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 62451262,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 5558695.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (sorted)",
            "value": 94310109,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/unsorted)",
            "value": 258489326.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 343487227,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 285978017,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/sorted)",
            "value": 120345497,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 274974480.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/set_elements!",
            "value": 2349708,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45392\nallocs=731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/foreach_pair (sorted)",
            "value": 69077658,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41696\nallocs=641\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 79317888.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/set_elements!",
            "value": 9141041,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45392\nallocs=731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/foreach_pair (sorted)",
            "value": 200767701,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41696\nallocs=641\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 232061496.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2837910861.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68954560\nallocs=2243875\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 28694713,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131888\nallocs=2540\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "66eba377d5c0f7f2463981f526f6027dd5477cc0",
          "message": "Docs: add logo",
          "timestamp": "2026-01-29T17:59:22+01:00",
          "tree_id": "b10d79d833831b3debbd5538654ab44802a753e7",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/66eba377d5c0f7f2463981f526f6027dd5477cc0"
        },
        "date": 1769706600068,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 13411829,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367520\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12759796,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 705892328.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5321312\nallocs=15977\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 588860245,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5214864\nallocs=15552\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8661451,
            "unit": "ns",
            "extra": "gctime=0\nmemory=213184\nallocs=2057\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1139543607,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9002808\nallocs=55594\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1382305443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42107504\nallocs=25046\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5520426.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73696\nallocs=1714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 873131,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84416\nallocs=2049\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2250595.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3220768439,
            "unit": "ns",
            "extra": "gctime=0\nmemory=334881968\nallocs=4179254\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 271806611,
            "unit": "ns",
            "extra": "gctime=27010534.5\nmemory=209976264\nallocs=7493721\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 7108406,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 1961163.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (sorted)",
            "value": 55255055,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/unsorted)",
            "value": 87508929,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 101669409,
            "unit": "ns",
            "extra": "gctime=0\nmemory=89602208\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 64336505,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/sorted)",
            "value": 75453277,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 64199372,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 6064163,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (sorted)",
            "value": 98958057,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/unsorted)",
            "value": 266665775,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 343019731,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 261514483,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/sorted)",
            "value": 107456195,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 260556080.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/set_elements!",
            "value": 2354244,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45392\nallocs=731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/foreach_pair (sorted)",
            "value": 63645847.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41696\nallocs=641\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 66958659,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/set_elements!",
            "value": 7980946,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45392\nallocs=731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/foreach_pair (sorted)",
            "value": 185704808,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41696\nallocs=641\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 215736848,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2781004656,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68889856\nallocs=2241403\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 28527685,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131888\nallocs=2540\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "juan-ignacio.polanco@cnrs.fr",
            "name": "Juan Ignacio Polanco",
            "username": "jipolanco"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": false,
          "id": "a6b650f10310c211848780db59b3b201936da80a",
          "message": "Update gif in README.md",
          "timestamp": "2026-01-30T12:10:26+01:00",
          "tree_id": "4c468662c667a27c4ace4092955a43bae935dff7",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/a6b650f10310c211848780db59b3b201936da80a"
        },
        "date": 1769772090352,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 13450451,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367520\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12799096.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 713552207,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5333408\nallocs=16181\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 596611194,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5220912\nallocs=15654\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8676443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=213184\nallocs=2057\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1170968768,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9057240\nallocs=56512\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1535702595.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42161936\nallocs=25964\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5501722,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73696\nallocs=1714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 872158,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84416\nallocs=2049\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2250579,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3235581719,
            "unit": "ns",
            "extra": "gctime=0\nmemory=334881968\nallocs=4179254\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 263581083.5,
            "unit": "ns",
            "extra": "gctime=26949325.5\nmemory=209976264\nallocs=7493721\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 6889657,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 1997469,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (sorted)",
            "value": 56286540,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/unsorted)",
            "value": 92218318,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 116481754,
            "unit": "ns",
            "extra": "gctime=0\nmemory=89602208\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 66094764,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/sorted)",
            "value": 78471339,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 65696439.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 6000130,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (sorted)",
            "value": 99332907.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/unsorted)",
            "value": 275330852,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 346399055,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 269226740,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/sorted)",
            "value": 109592187.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 268536967,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/set_elements!",
            "value": 2345357.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45392\nallocs=731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/foreach_pair (sorted)",
            "value": 63720501,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41696\nallocs=641\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 65387712,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/set_elements!",
            "value": 7664061,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45392\nallocs=731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/foreach_pair (sorted)",
            "value": 187939010,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41696\nallocs=641\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 214747624.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2747229380.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68920096\nallocs=2241913\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 27783961,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131888\nallocs=2540\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "5dd7a432941627342c34e288346389db7040ea9f",
          "message": "Docs: add note on heterogeneous clusters",
          "timestamp": "2026-01-30T13:11:07+01:00",
          "tree_id": "7101c56677ae4f0d31e581537d4ff9147aafa897",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/5dd7a432941627342c34e288346389db7040ea9f"
        },
        "date": 1769775720058,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 13503112,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367520\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12937791,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 725800581,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5321312\nallocs=15977\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 606797426,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5214864\nallocs=15552\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8702284,
            "unit": "ns",
            "extra": "gctime=0\nmemory=213184\nallocs=2057\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1177986419,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9002808\nallocs=55594\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1629628455.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42107504\nallocs=25046\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5557684,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73696\nallocs=1714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 875261,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84416\nallocs=2049\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2254820,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3219958334.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=334883072\nallocs=4179267\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 272609774,
            "unit": "ns",
            "extra": "gctime=27894789\nmemory=209976264\nallocs=7493721\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 7179932,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 2082059,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (sorted)",
            "value": 55918147,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/unsorted)",
            "value": 88429897.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 108857893,
            "unit": "ns",
            "extra": "gctime=0\nmemory=89602208\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 68292124,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/sorted)",
            "value": 644106431,
            "unit": "ns",
            "extra": "gctime=290934159\nmemory=1568660416\nallocs=14655084\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 80904633,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 6090234,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (sorted)",
            "value": 98630421,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/unsorted)",
            "value": 276066567.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 365631657,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 277968180,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/sorted)",
            "value": 110748892,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 267490270,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/set_elements!",
            "value": 2362373.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45392\nallocs=731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/foreach_pair (sorted)",
            "value": 62694007.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41696\nallocs=641\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 64052860,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/set_elements!",
            "value": 7651708,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45392\nallocs=731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/foreach_pair (sorted)",
            "value": 185920693,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41696\nallocs=641\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 220349226,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2757033208.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68884704\nallocs=2241091\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 26484448,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131888\nallocs=2540\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "fe398de991f97236759a6f9e5b0defff9c352ab5",
          "message": "Docs: add favicon",
          "timestamp": "2026-01-30T14:36:54+01:00",
          "tree_id": "83815b925af23f2f4873e928e0b803a42b5e1369",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/fe398de991f97236759a6f9e5b0defff9c352ab5"
        },
        "date": 1769780848292,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 13287498,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367520\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12755307,
            "unit": "ns",
            "extra": "gctime=0\nmemory=726528\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 722719070,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5330752\nallocs=16279\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 596833142,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5218256\nallocs=15752\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8673658,
            "unit": "ns",
            "extra": "gctime=0\nmemory=213184\nallocs=2057\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1157757267,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9057240\nallocs=56512\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1542177305,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42161936\nallocs=25964\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5480969,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73696\nallocs=1714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 872129,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84416\nallocs=2049\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2308037,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3276663226,
            "unit": "ns",
            "extra": "gctime=0\nmemory=334882224\nallocs=4179258\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 257049328,
            "unit": "ns",
            "extra": "gctime=26572640\nmemory=209976264\nallocs=7493721\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 6971057.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 2013386,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (sorted)",
            "value": 56687794,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/unsorted)",
            "value": 87781221,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 108935314,
            "unit": "ns",
            "extra": "gctime=0\nmemory=89602208\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 64759735,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/sorted)",
            "value": 76945856,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 64320117,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 5689487.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (sorted)",
            "value": 96022257,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/unsorted)",
            "value": 254155781,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 335561114,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 257650675,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/sorted)",
            "value": 104802529.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 254796723.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/set_elements!",
            "value": 2325876.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45392\nallocs=731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/foreach_pair (sorted)",
            "value": 62978340,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41696\nallocs=641\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 63385880,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/set_elements!",
            "value": 7246826.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45392\nallocs=731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/foreach_pair (sorted)",
            "value": 184626855,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41696\nallocs=641\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 209320959.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2741477197.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68909328\nallocs=2242304\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 27087060,
            "unit": "ns",
            "extra": "gctime=0\nmemory=131888\nallocs=2540\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "juan-ignacio.polanco@cnrs.fr",
            "name": "Juan Ignacio Polanco",
            "username": "jipolanco"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "42f51429b56d707be0fb30e0962956f266f4172d",
          "message": "Implement HostVector for more robust CPU <-> GPU data transfers (#93)\n\n* Define copy_device_to_host!\n\n* Implement HostVector + package extensions\n\n* Use HostVector in host-device copies (untested)\n\n* Update Project.toml\n\n* Fix HostVector constructor\n\n* Fix resize_no_copy!\n\n* Make things work\n\n* Possible improvements\n\n* Fix aliasing issue\n\n* Make things work on CUDA _and_ OpenCL\n\n* Add some @debug output\n\n* ring_collision test: avoid scalar indexing on GPU",
          "timestamp": "2026-02-03T10:54:04+01:00",
          "tree_id": "0efc68e64ab2c7a0fd9dd2346f5f52f2a98eafb7",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/42f51429b56d707be0fb30e0962956f266f4172d"
        },
        "date": 1770113528390,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 12987676,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367264\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12099436,
            "unit": "ns",
            "extra": "gctime=0\nmemory=710288\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 747874857,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5234768\nallocs=14573\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 610072790,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5128640\nallocs=14148\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8373264,
            "unit": "ns",
            "extra": "gctime=0\nmemory=213056\nallocs=2049\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1297940588.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9000504\nallocs=55450\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1588863686.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42107504\nallocs=25046\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5371297.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73696\nallocs=1714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 904724,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84416\nallocs=2049\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2391817,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 2911450515.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=334882224\nallocs=4179258\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 297556853,
            "unit": "ns",
            "extra": "gctime=25337870\nmemory=209976264\nallocs=7493721\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 7484561.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 3703703,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (sorted)",
            "value": 63251403,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/unsorted)",
            "value": 103825293.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 96848354,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 79583698,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/sorted)",
            "value": 514688732,
            "unit": "ns",
            "extra": "gctime=220783885.5\nmemory=1341608416\nallocs=10481334\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 79777153,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 6538235,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (sorted)",
            "value": 106398912,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/unsorted)",
            "value": 187161838,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 270675158,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 172246024.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/sorted)",
            "value": 119602623,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 172025749,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/set_elements!",
            "value": 3031551,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45392\nallocs=731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/foreach_pair (sorted)",
            "value": 79075316,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41696\nallocs=641\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 76999000,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/set_elements!",
            "value": 8842278,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45392\nallocs=731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/foreach_pair (sorted)",
            "value": 250311432,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41696\nallocs=641\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 177371558,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2908501792,
            "unit": "ns",
            "extra": "gctime=0\nmemory=69458712\nallocs=2236457\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 25584792,
            "unit": "ns",
            "extra": "gctime=0\nmemory=105936\nallocs=2011\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "e27b3a4caf12b430283083a61da3dd1f1d814b64",
          "message": "Update CHANGELOG",
          "timestamp": "2026-02-03T11:45:03+01:00",
          "tree_id": "0901f06d1d876a78def067bfcb79fe3773f57f6a",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/e27b3a4caf12b430283083a61da3dd1f1d814b64"
        },
        "date": 1770116151918,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 13186510,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367264\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12247388,
            "unit": "ns",
            "extra": "gctime=0\nmemory=710288\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 766010613,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5252912\nallocs=14879\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 626182009,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5140736\nallocs=14352\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8727695,
            "unit": "ns",
            "extra": "gctime=0\nmemory=213184\nallocs=2057\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1336464719.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9057240\nallocs=56512\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1599909352.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42161936\nallocs=25964\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5404500,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73696\nallocs=1714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 896768.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84416\nallocs=2049\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2452036,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3135849435,
            "unit": "ns",
            "extra": "gctime=0\nmemory=334882480\nallocs=4179262\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 301614223,
            "unit": "ns",
            "extra": "gctime=25544467\nmemory=209976264\nallocs=7493721\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 7590537.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 3801886,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (sorted)",
            "value": 66964370,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/unsorted)",
            "value": 110620869,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 101086217,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 83382657.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/sorted)",
            "value": 91222103,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 81829340,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 10430149,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (sorted)",
            "value": 117846254,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/unsorted)",
            "value": 250631000.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 349519779,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 210634193,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/sorted)",
            "value": 126172199,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 211115097,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/set_elements!",
            "value": 3252830,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45392\nallocs=731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/foreach_pair (sorted)",
            "value": 83609694.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41696\nallocs=641\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 83345745,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/set_elements!",
            "value": 11580553,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45392\nallocs=731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/foreach_pair (sorted)",
            "value": 259529569,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41696\nallocs=641\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 223529360,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2935463258,
            "unit": "ns",
            "extra": "gctime=0\nmemory=69414264\nallocs=2234903\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 26006023,
            "unit": "ns",
            "extra": "gctime=0\nmemory=105936\nallocs=2011\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "juan-ignacio.polanco@cnrs.fr",
            "name": "Juan Ignacio Polanco",
            "username": "jipolanco"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "7e8beb943c14761ad4a4a38490a35a51766556ea",
          "message": "Parallelise some CPU -> CPU copies (#94)\n\n* Parallelise some CPU -> CPU copies\n\n* Fix GPU->CPU transfer",
          "timestamp": "2026-02-03T13:50:39+01:00",
          "tree_id": "ad2847287abb75a48fcaf17a292d5242d45ef088",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/7e8beb943c14761ad4a4a38490a35a51766556ea"
        },
        "date": 1770124172604,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 13341600,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367264\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12337965,
            "unit": "ns",
            "extra": "gctime=0\nmemory=710288\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 788101816,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5286960\nallocs=15303\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 644903313,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5175904\nallocs=14810\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8793953,
            "unit": "ns",
            "extra": "gctime=0\nmemory=213184\nallocs=2057\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1410630524,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9024984\nallocs=55900\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1827711803,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42129680\nallocs=25352\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5459207,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4016\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 902573,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84416\nallocs=2049\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2624899.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3046934349,
            "unit": "ns",
            "extra": "gctime=0\nmemory=334882736\nallocs=4179266\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 310750906,
            "unit": "ns",
            "extra": "gctime=21504227.5\nmemory=209976264\nallocs=7493721\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 7942486,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 4571022,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (sorted)",
            "value": 74274515.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/unsorted)",
            "value": 140387535,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 124568470,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 96400011,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/sorted)",
            "value": 633341819.5,
            "unit": "ns",
            "extra": "gctime=262953522.5\nmemory=1568662320\nallocs=14655119\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 98167537,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 11918685,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (sorted)",
            "value": 122906927,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/unsorted)",
            "value": 265748452,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 377938022,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 243470786,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/sorted)",
            "value": 140005118.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 232929710.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/set_elements!",
            "value": 3872259,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45392\nallocs=731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/foreach_pair (sorted)",
            "value": 88902224,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41696\nallocs=641\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 88637624,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/set_elements!",
            "value": 12396516.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45392\nallocs=731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/foreach_pair (sorted)",
            "value": 262328665,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41696\nallocs=641\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 237055074,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 3047789015.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68753440\nallocs=2237451\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 26897078,
            "unit": "ns",
            "extra": "gctime=0\nmemory=124688\nallocs=2277\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "dbbd29218d5d247bbb53e6266b9f5b36a4f04d11",
          "message": "v0.32.18",
          "timestamp": "2026-02-03T13:54:38+01:00",
          "tree_id": "c090a566eac2ae4f6027f188187fdd30bfe82a84",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/dbbd29218d5d247bbb53e6266b9f5b36a4f04d11"
        },
        "date": 1770124419545,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 13705909.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4367264\nallocs=7038\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 13027046,
            "unit": "ns",
            "extra": "gctime=0\nmemory=710288\nallocs=6035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 715094549,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5286960\nallocs=15303\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 594613633,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5175904\nallocs=14810\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8901223.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=213184\nallocs=2057\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1150799312,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9024984\nallocs=55900\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1614336174,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42129680\nallocs=25352\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5510932,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73696\nallocs=1714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 878925,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84416\nallocs=2049\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2262900,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3288323011,
            "unit": "ns",
            "extra": "gctime=0\nmemory=334881456\nallocs=4179246\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 273811450.5,
            "unit": "ns",
            "extra": "gctime=26365655\nmemory=209976264\nallocs=7493721\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 7127167,
            "unit": "ns",
            "extra": "gctime=0\nmemory=275376\nallocs=5593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 1981180,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (sorted)",
            "value": 61467128,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/unsorted)",
            "value": 115982927,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 137459981,
            "unit": "ns",
            "extra": "gctime=0\nmemory=89602208\nallocs=1200022\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 87479943,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/sorted)",
            "value": 85321061,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 84549575,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 7762023.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (sorted)",
            "value": 112598810,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/unsorted)",
            "value": 317517549,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 402924267,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 306936114,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/sorted)",
            "value": 123446293,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 312402902,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/set_elements!",
            "value": 2369327,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45392\nallocs=731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/foreach_pair (sorted)",
            "value": 70643502,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41696\nallocs=641\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 87415538,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/set_elements!",
            "value": 8941083,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45392\nallocs=731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/foreach_pair (sorted)",
            "value": 195687269,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41696\nallocs=641\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/OpenCLBackend/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 250570997.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9328\nallocs=122\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2775158505.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=69990392\nallocs=2237519\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 22937974,
            "unit": "ns",
            "extra": "gctime=0\nmemory=124384\nallocs=2271\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "814e80f613cf1ceacc2a010c58fc259cb3d432fa",
          "message": "Disable OpenCL benchmarks",
          "timestamp": "2026-02-24T17:25:50+01:00",
          "tree_id": "c4509f6689f5e484fb0632c568523c094dde3ed7",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/814e80f613cf1ceacc2a010c58fc259cb3d432fa"
        },
        "date": 1771951353254,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 12257811,
            "unit": "ns",
            "extra": "gctime=0\nmemory=18528\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12723522,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4448\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 702517002,
            "unit": "ns",
            "extra": "gctime=0\nmemory=232384\nallocs=2298\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 586220131,
            "unit": "ns",
            "extra": "gctime=0\nmemory=185520\nallocs=1805\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8635731,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4016\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1145134783,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4681080\nallocs=7516\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1470889322.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42129680\nallocs=25352\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5495785,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4016\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 827553.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4016\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2259717.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3180674512.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68636176\nallocs=176606\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 208077289,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9931848\nallocs=2697\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 6817192,
            "unit": "ns",
            "extra": "gctime=0\nmemory=99376\nallocs=2593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 1999610.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (sorted)",
            "value": 57545147,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7136\nallocs=74\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/unsorted)",
            "value": 96663299.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 80270925,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 70892452.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/sorted)",
            "value": 81845768.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7136\nallocs=74\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 94783942,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 6935452,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (sorted)",
            "value": 105121312,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7136\nallocs=74\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/unsorted)",
            "value": 289545564,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 371330530.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 285712584.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/sorted)",
            "value": 113796532,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7136\nallocs=74\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 291280078.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2714639836.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4421832\nallocs=10963\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 22907973,
            "unit": "ns",
            "extra": "gctime=0\nmemory=29248\nallocs=394\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "distinct": false,
          "id": "814e80f613cf1ceacc2a010c58fc259cb3d432fa",
          "message": "Disable OpenCL benchmarks",
          "timestamp": "2026-02-24T17:25:50+01:00",
          "tree_id": "c4509f6689f5e484fb0632c568523c094dde3ed7",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/814e80f613cf1ceacc2a010c58fc259cb3d432fa"
        },
        "date": 1771952118698,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 12207435,
            "unit": "ns",
            "extra": "gctime=0\nmemory=18528\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12730851,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4448\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 702072488.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=232384\nallocs=2298\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 604707396,
            "unit": "ns",
            "extra": "gctime=0\nmemory=185520\nallocs=1805\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8882569.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4016\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1155030698,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4681080\nallocs=7516\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1547752785,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42129680\nallocs=25352\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5491986,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4016\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 877647,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4016\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2267466,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3143765108.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68636176\nallocs=176606\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 200383066,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9931848\nallocs=2697\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 6816905,
            "unit": "ns",
            "extra": "gctime=0\nmemory=99376\nallocs=2593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 1939464,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (sorted)",
            "value": 61372180,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/unsorted)",
            "value": 121494297,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 95258875,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 82379510,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/sorted)",
            "value": 95879812.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 92580244,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 8508091,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (sorted)",
            "value": 123746389.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/unsorted)",
            "value": 344600491,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 439313618.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 332917068,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/sorted)",
            "value": 143726567.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 308140112,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2715194554.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4421832\nallocs=10963\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 22941227,
            "unit": "ns",
            "extra": "gctime=0\nmemory=29248\nallocs=394\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "juan-ignacio.polanco@cnrs.fr",
            "name": "Juan Ignacio Polanco",
            "username": "jipolanco"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "cf5f8f3784cc60ad693faff55e6e69aafcb6780e",
          "message": "Implement Ewald summation based on Kaiser-Bessel kernel (#95)\n\n* Define and use GaussianMollifier\n\n* shortrange -> longrange\n\n* Fix loading\n\n* Fix SIMD version\n\n* Allow α = Zero() in GaussianMollifier\n\n* Rename mollifier -> splitting\n\n* Replace α with GaussianSplitting in params\n\n* Update some docs\n\n* Rename files: \"mollifiers\" -> \"splitting\"\n\n* Allow splitting argument in ParamsBiotSavart\n\n* Fix GaussianSplitting constructor\n\n* Use splitting argument in tutorials\n\n* Update docs\n\n* Add more detailed docs for GaussianSplitting\n\n* Start work on KB splitting\n\n* Update docs\n\n* Add Chebyshev approximation of KB kernels\n\n* Add ChebyshevSeries tests\n\n* Allow KB on GPU + estimate background constant\n\n* Fix long-range correction\n\n* Update trefoil test\n\n* Update some tests\n\n* Test and fix background vorticity correction\n\n* Fix vortex ring tutorial\n\n* Update leapfrogging test\n\n* Remove BiotSavart.autotune\n\n* Test KB with OpenCLBackend\n\n* Fix forcing test\n\n* Use KaiserBesselSplitting in forcing tests\n\n* Add docs/package.json\n\n* vitepress: update config.mts\n\n* Update from vitepress template\n\n* Docs: use mathjax-pagella-font\n\n* Update benchmarks\n\n* Truncate Chebyshev series earlier for β < 19\n\n* Update docs\n\n* Update KW tests + precompilation",
          "timestamp": "2026-02-25T12:21:52+01:00",
          "tree_id": "91cf8f3fc95795e091116617f6aea3dd025a9cf8",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/cf5f8f3784cc60ad693faff55e6e69aafcb6780e"
        },
        "date": 1772019810339,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 13714966,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20320\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting{Float64, 3}/add_local_integrals!",
            "value": 13772684.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20320\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting{Float64, 3}/velocity + streamfunction",
            "value": 847905246.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=238944\nallocs=2498\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting{Float64, 3}/velocity",
            "value": 691439787.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=186016\nallocs=1903\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting{Float64, 3, VortexPasta.BiotSavart.ChebyshevApproximations.ChebyshevSeries{:even, Float64, Vector{Float64}}, VortexPasta.BiotSavart.ChebyshevApproximations.ChebyshevSeries{:odd, Float64, Vector{Float64}}}/add_local_integrals!",
            "value": 13715683,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22368\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting{Float64, 3, VortexPasta.BiotSavart.ChebyshevApproximations.ChebyshevSeries{:even, Float64, Vector{Float64}}, VortexPasta.BiotSavart.ChebyshevApproximations.ChebyshevSeries{:odd, Float64, Vector{Float64}}}/velocity + streamfunction",
            "value": 645014314,
            "unit": "ns",
            "extra": "gctime=0\nmemory=242640\nallocs=2498\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting{Float64, 3, VortexPasta.BiotSavart.ChebyshevApproximations.ChebyshevSeries{:even, Float64, Vector{Float64}}, VortexPasta.BiotSavart.ChebyshevApproximations.ChebyshevSeries{:odd, Float64, Vector{Float64}}}/velocity",
            "value": 527151706.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=189536\nallocs=1903\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 13985499,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5088\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 823358032,
            "unit": "ns",
            "extra": "gctime=0\nmemory=238944\nallocs=2498\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 697437733,
            "unit": "ns",
            "extra": "gctime=0\nmemory=186016\nallocs=1903\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 9659362,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1272854948,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4744728\nallocs=8434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1609418422,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42184688\nallocs=26270\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 6165905,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 934107.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2372924,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3497538977,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68636688\nallocs=176606\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 182568499,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9931848\nallocs=2697\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 7220383,
            "unit": "ns",
            "extra": "gctime=0\nmemory=99376\nallocs=2593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 2042200,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (sorted)",
            "value": 65944870.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/unsorted)",
            "value": 97952927.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 91329053.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 86602507,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/sorted)",
            "value": 87262322,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 78004507,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 7522697,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (sorted)",
            "value": 123691463.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/unsorted)",
            "value": 263955098,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 442797269.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 349154528,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/sorted)",
            "value": 126590941,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 358697121.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 3125021900,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4164952\nallocs=11500\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 24703947,
            "unit": "ns",
            "extra": "gctime=0\nmemory=29248\nallocs=394\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "9226d8d74fae8c9c84f1e29ff1fd164147893cb3",
          "message": "Try to fix precompilation on Julia 1.11",
          "timestamp": "2026-02-25T12:45:03+01:00",
          "tree_id": "d70a6640fe49c2c58aa6e7a829c7fac9b427a937",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/9226d8d74fae8c9c84f1e29ff1fd164147893cb3"
        },
        "date": 1772021163509,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 12398630,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20320\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting{Float64, 3}/add_local_integrals!",
            "value": 12403150.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20320\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting{Float64, 3}/velocity + streamfunction",
            "value": 756525600,
            "unit": "ns",
            "extra": "gctime=0\nmemory=232896\nallocs=2396\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting{Float64, 3}/velocity",
            "value": 621399964,
            "unit": "ns",
            "extra": "gctime=0\nmemory=186016\nallocs=1903\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting{Float64, 3, VortexPasta.BiotSavart.ChebyshevApproximations.ChebyshevSeries{:even, Float64, Vector{Float64}}, VortexPasta.BiotSavart.ChebyshevApproximations.ChebyshevSeries{:odd, Float64, Vector{Float64}}}/add_local_integrals!",
            "value": 14392666.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22368\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting{Float64, 3, VortexPasta.BiotSavart.ChebyshevApproximations.ChebyshevSeries{:even, Float64, Vector{Float64}}, VortexPasta.BiotSavart.ChebyshevApproximations.ChebyshevSeries{:odd, Float64, Vector{Float64}}}/velocity + streamfunction",
            "value": 628460509,
            "unit": "ns",
            "extra": "gctime=0\nmemory=236592\nallocs=2396\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting{Float64, 3, VortexPasta.BiotSavart.ChebyshevApproximations.ChebyshevSeries{:even, Float64, Vector{Float64}}, VortexPasta.BiotSavart.ChebyshevApproximations.ChebyshevSeries{:odd, Float64, Vector{Float64}}}/velocity",
            "value": 495827104,
            "unit": "ns",
            "extra": "gctime=0\nmemory=189536\nallocs=1903\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12119918.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5088\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 762731206,
            "unit": "ns",
            "extra": "gctime=0\nmemory=232896\nallocs=2396\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 622863860,
            "unit": "ns",
            "extra": "gctime=0\nmemory=186016\nallocs=1903\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8494564,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1331885199,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4690296\nallocs=7516\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1750335479,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42130256\nallocs=25352\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5355066,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 822391.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2346226,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3172251508.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68636176\nallocs=176606\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 235503102,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9931848\nallocs=2697\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 7591485.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=99376\nallocs=2593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 3629262,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (sorted)",
            "value": 67754082.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/unsorted)",
            "value": 106162023,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 102436721,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 78473130.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/sorted)",
            "value": 109425826,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 96398867,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 10260628,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (sorted)",
            "value": 112709164.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/unsorted)",
            "value": 212748592,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 308317036,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 187067599,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/sorted)",
            "value": 121733982,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 193032359,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2904257748.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4158904\nallocs=11398\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 26156413,
            "unit": "ns",
            "extra": "gctime=0\nmemory=29248\nallocs=394\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "fff87316e3b9b99da813d76f9097e7ebbad227ff",
          "message": "Benchmarks: fix names of splitting types",
          "timestamp": "2026-02-25T13:15:59+01:00",
          "tree_id": "5e914461a11c9c5df55b8593eda4608f73e529b2",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/fff87316e3b9b99da813d76f9097e7ebbad227ff"
        },
        "date": 1772022837084,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/KaiserBesselSplitting/add_local_integrals!",
            "value": 12232626,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22368\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting/velocity + streamfunction",
            "value": 562562975,
            "unit": "ns",
            "extra": "gctime=0\nmemory=230544\nallocs=2294\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting/velocity",
            "value": 440810948.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=183488\nallocs=1801\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 12272082,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20320\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12635361.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5088\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 705020839,
            "unit": "ns",
            "extra": "gctime=0\nmemory=226848\nallocs=2294\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 585257357,
            "unit": "ns",
            "extra": "gctime=0\nmemory=179968\nallocs=1801\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/add_local_integrals!",
            "value": 12267004.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20320\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/velocity + streamfunction",
            "value": 710997800,
            "unit": "ns",
            "extra": "gctime=0\nmemory=226848\nallocs=2294\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/velocity",
            "value": 582146400,
            "unit": "ns",
            "extra": "gctime=0\nmemory=179968\nallocs=1801\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8775133.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1178981039,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4690296\nallocs=7516\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1660216304,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42130256\nallocs=25352\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5465863,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 880055,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2221733,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3157543482,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68636432\nallocs=176606\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 216401936,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9931848\nallocs=2697\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 7043230.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=99376\nallocs=2593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 1967490,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (sorted)",
            "value": 59389513,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/unsorted)",
            "value": 96860190.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 81073597,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 65962437.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/sorted)",
            "value": 78189171,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 64039452.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 6675229.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (sorted)",
            "value": 103029249,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/unsorted)",
            "value": 300332422,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 374978055,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 294877297.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/sorted)",
            "value": 114967808,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 286802246,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2729948476,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4405144\nallocs=10947\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 22886198,
            "unit": "ns",
            "extra": "gctime=0\nmemory=29248\nallocs=394\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "55aad2da7a046c0d901ccdac3decf64a7527a3dc",
          "message": "v0.33.0",
          "timestamp": "2026-02-25T13:37:39+01:00",
          "tree_id": "43b641f1dbd41b25b6f0a83c5fa2a770f1c41aea",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/55aad2da7a046c0d901ccdac3decf64a7527a3dc"
        },
        "date": 1772024099183,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/KaiserBesselSplitting/add_local_integrals!",
            "value": 12170878,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22368\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting/velocity + streamfunction",
            "value": 617793977,
            "unit": "ns",
            "extra": "gctime=0\nmemory=242640\nallocs=2498\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting/velocity",
            "value": 484696925,
            "unit": "ns",
            "extra": "gctime=0\nmemory=189536\nallocs=1903\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 12394470,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20320\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 11899805,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5088\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 755198961,
            "unit": "ns",
            "extra": "gctime=0\nmemory=238944\nallocs=2498\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 619405842,
            "unit": "ns",
            "extra": "gctime=0\nmemory=186016\nallocs=1903\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/add_local_integrals!",
            "value": 12157746,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20320\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/velocity + streamfunction",
            "value": 750187641,
            "unit": "ns",
            "extra": "gctime=0\nmemory=238944\nallocs=2498\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/velocity",
            "value": 619926867,
            "unit": "ns",
            "extra": "gctime=0\nmemory=186016\nallocs=1903\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8371878,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1338793348.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4744728\nallocs=8434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1598303776.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42184688\nallocs=26270\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5320776.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 890967,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2417391,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3315873470,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68636432\nallocs=176606\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 232098624.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9931848\nallocs=2697\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 7406620.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=99376\nallocs=2593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 3587808.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (sorted)",
            "value": 61984144,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/unsorted)",
            "value": 106152070,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 97532359,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 80013158,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/sorted)",
            "value": 85142987.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 80203944,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 6172579,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (sorted)",
            "value": 106384625,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/unsorted)",
            "value": 190028877,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 262842672,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 170473683,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/sorted)",
            "value": 121896501,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 169826184,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2944661462.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4164840\nallocs=11499\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 25840118,
            "unit": "ns",
            "extra": "gctime=0\nmemory=29248\nallocs=394\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "7a85f098ff74b70b2302e380a3543fbf0ba8e4b5",
          "message": "v0.33.1",
          "timestamp": "2026-02-25T15:58:45+01:00",
          "tree_id": "f6f6ae16e0921b7fce66a35f7db0519853c67ea4",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/7a85f098ff74b70b2302e380a3543fbf0ba8e4b5"
        },
        "date": 1772032568441,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/KaiserBesselSplitting/add_local_integrals!",
            "value": 12218221.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22368\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting/velocity + streamfunction",
            "value": 538727539.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=230544\nallocs=2294\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting/velocity",
            "value": 427775216.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=183488\nallocs=1801\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 12197119.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20320\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12519507,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5088\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 691760374,
            "unit": "ns",
            "extra": "gctime=0\nmemory=226848\nallocs=2294\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 571776742,
            "unit": "ns",
            "extra": "gctime=0\nmemory=179968\nallocs=1801\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/add_local_integrals!",
            "value": 12166342.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20320\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/velocity + streamfunction",
            "value": 691341839.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=226848\nallocs=2294\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/velocity",
            "value": 570827914,
            "unit": "ns",
            "extra": "gctime=0\nmemory=179968\nallocs=1801\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8627772.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1136408781,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4690296\nallocs=7516\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1413780510,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42130256\nallocs=25352\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5481889,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 824910,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2219735,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3154751636,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68636176\nallocs=176606\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 213757498.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9931848\nallocs=2697\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 6749345,
            "unit": "ns",
            "extra": "gctime=0\nmemory=99376\nallocs=2593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 1990360,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (sorted)",
            "value": 53759594,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/unsorted)",
            "value": 88570388,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 76232482.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 63925597,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/sorted)",
            "value": 74033452,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 63556789,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 5452720,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (sorted)",
            "value": 95607099.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/unsorted)",
            "value": 255532547.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 325982920,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 251886230,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/sorted)",
            "value": 108474446,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 252975580.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2718462346.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4405144\nallocs=10947\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 22001164,
            "unit": "ns",
            "extra": "gctime=0\nmemory=29248\nallocs=394\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "49699333+dependabot[bot]@users.noreply.github.com",
            "name": "dependabot[bot]",
            "username": "dependabot[bot]"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "5605689c1ae4829c15598bbe933b3fce560b9f7e",
          "message": "Bump actions/upload-artifact from 6 to 7 (#96)\n\nBumps [actions/upload-artifact](https://github.com/actions/upload-artifact) from 6 to 7.\n- [Release notes](https://github.com/actions/upload-artifact/releases)\n- [Commits](https://github.com/actions/upload-artifact/compare/v6...v7)\n\n---\nupdated-dependencies:\n- dependency-name: actions/upload-artifact\n  dependency-version: '7'\n  dependency-type: direct:production\n  update-type: version-update:semver-major\n...\n\nSigned-off-by: dependabot[bot] <support@github.com>\nCo-authored-by: dependabot[bot] <49699333+dependabot[bot]@users.noreply.github.com>",
          "timestamp": "2026-03-02T09:07:08+01:00",
          "tree_id": "5d41d3e09ee6d37a403b6f9dff126fca8df4c84a",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/5605689c1ae4829c15598bbe933b3fce560b9f7e"
        },
        "date": 1772439428451,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/KaiserBesselSplitting/add_local_integrals!",
            "value": 12286571.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22368\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting/velocity + streamfunction",
            "value": 562011969,
            "unit": "ns",
            "extra": "gctime=0\nmemory=236592\nallocs=2396\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting/velocity",
            "value": 435400941.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=189536\nallocs=1903\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 12332146,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20320\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12754853.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5088\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 720908384,
            "unit": "ns",
            "extra": "gctime=0\nmemory=232896\nallocs=2396\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 584779487,
            "unit": "ns",
            "extra": "gctime=0\nmemory=186016\nallocs=1903\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/add_local_integrals!",
            "value": 12268870,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20320\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/velocity + streamfunction",
            "value": 713773833,
            "unit": "ns",
            "extra": "gctime=0\nmemory=232896\nallocs=2396\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/velocity",
            "value": 587969426,
            "unit": "ns",
            "extra": "gctime=0\nmemory=186016\nallocs=1903\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8688928.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1160735845,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4690296\nallocs=7516\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1476293676.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42130256\nallocs=25352\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5498806,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 812881,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2221435,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 4029128025.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68636176\nallocs=176606\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 229044907.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9931848\nallocs=2697\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 7110304,
            "unit": "ns",
            "extra": "gctime=0\nmemory=99376\nallocs=2593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 2083166.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (sorted)",
            "value": 59296177,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/unsorted)",
            "value": 90017249,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 79172393.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 72931654,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/sorted)",
            "value": 84977743,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 71609809.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 6141684,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (sorted)",
            "value": 94567962,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/unsorted)",
            "value": 266591558,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 373765308,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 260670958,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/sorted)",
            "value": 111100572,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 264341923,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2746502907.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4158792\nallocs=11397\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 23846756,
            "unit": "ns",
            "extra": "gctime=0\nmemory=29248\nallocs=394\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "49699333+dependabot[bot]@users.noreply.github.com",
            "name": "dependabot[bot]",
            "username": "dependabot[bot]"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "89a4f6febb0be9efa565bb685be1ff6b66b6e47e",
          "message": "Bump julia-actions/cache from 2 to 3 (#97)\n\nBumps [julia-actions/cache](https://github.com/julia-actions/cache) from 2 to 3.\n- [Release notes](https://github.com/julia-actions/cache/releases)\n- [Commits](https://github.com/julia-actions/cache/compare/v2...v3)\n\n---\nupdated-dependencies:\n- dependency-name: julia-actions/cache\n  dependency-version: '3'\n  dependency-type: direct:production\n  update-type: version-update:semver-major\n...\n\nSigned-off-by: dependabot[bot] <support@github.com>\nCo-authored-by: dependabot[bot] <49699333+dependabot[bot]@users.noreply.github.com>",
          "timestamp": "2026-03-09T09:47:30+01:00",
          "tree_id": "b9fe28c07fea72b4919bc135707fcdf8110c0c70",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/89a4f6febb0be9efa565bb685be1ff6b66b6e47e"
        },
        "date": 1773047146915,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/KaiserBesselSplitting/add_local_integrals!",
            "value": 12478950,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22368\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting/velocity + streamfunction",
            "value": 676305803,
            "unit": "ns",
            "extra": "gctime=0\nmemory=248688\nallocs=2600\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting/velocity",
            "value": 537410087,
            "unit": "ns",
            "extra": "gctime=0\nmemory=195584\nallocs=2005\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 12447483,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20320\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12138192,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5088\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 793966012,
            "unit": "ns",
            "extra": "gctime=0\nmemory=244992\nallocs=2600\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 651475499,
            "unit": "ns",
            "extra": "gctime=0\nmemory=192064\nallocs=2005\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/add_local_integrals!",
            "value": 12330193,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20320\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/velocity + streamfunction",
            "value": 825856994,
            "unit": "ns",
            "extra": "gctime=0\nmemory=244992\nallocs=2600\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/velocity",
            "value": 669146027.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=192064\nallocs=2005\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8510383.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1383459139,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4744728\nallocs=8434\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1815565332,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42184688\nallocs=26270\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5357104,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 880780,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2455846,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3951848667.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68636432\nallocs=176606\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 244477122,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9931848\nallocs=2697\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 8283909,
            "unit": "ns",
            "extra": "gctime=0\nmemory=99376\nallocs=2593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 4145753,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (sorted)",
            "value": 71422595.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/unsorted)",
            "value": 119195100,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 117085147,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 90077942,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/sorted)",
            "value": 96599712,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 93948491,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 10977410,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (sorted)",
            "value": 125329521.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/unsorted)",
            "value": 276351622,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 365714450,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 222410174,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/sorted)",
            "value": 142898439,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 230717393,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 3107236652.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3318064\nallocs=12073\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 28899150.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=29248\nallocs=394\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "9903a04e2d0a62b83dc07f65451022590eca9474",
          "message": "Make summary(::ParamsBiotSavart) more compact\n\nIn particular, this improves readability of `print(::VortexFilamentProblem)`.",
          "timestamp": "2026-03-19T12:55:15+01:00",
          "tree_id": "1020b9e09d055cd02691f3c3949f340726b683ce",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/9903a04e2d0a62b83dc07f65451022590eca9474"
        },
        "date": 1773922401033,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/KaiserBesselSplitting/add_local_integrals!",
            "value": 12141063,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22368\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting/velocity + streamfunction",
            "value": 562976564,
            "unit": "ns",
            "extra": "gctime=0\nmemory=230544\nallocs=2294\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting/velocity",
            "value": 445576390.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=183488\nallocs=1801\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 12225674,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20320\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12529759.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5088\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 717386922,
            "unit": "ns",
            "extra": "gctime=0\nmemory=226848\nallocs=2294\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 591246498,
            "unit": "ns",
            "extra": "gctime=0\nmemory=179968\nallocs=1801\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/add_local_integrals!",
            "value": 12230172.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20320\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/velocity + streamfunction",
            "value": 718973501,
            "unit": "ns",
            "extra": "gctime=0\nmemory=226848\nallocs=2294\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/velocity",
            "value": 601849678,
            "unit": "ns",
            "extra": "gctime=0\nmemory=179968\nallocs=1801\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8640932,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1177878198,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4690296\nallocs=7516\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1457658511.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42130256\nallocs=25352\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5491223,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 818045,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2230254,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 4065421151.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68636176\nallocs=176606\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 225160548,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9931848\nallocs=2697\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 7062289,
            "unit": "ns",
            "extra": "gctime=0\nmemory=99376\nallocs=2593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 2143049.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (sorted)",
            "value": 58927396,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/unsorted)",
            "value": 101704834,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 89911202.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 72999547,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/sorted)",
            "value": 83230578.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 73214525,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 6919631,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (sorted)",
            "value": 101676407,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/unsorted)",
            "value": 281206788,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 332704568,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 256848275,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/sorted)",
            "value": 110641908.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 273519427,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2764522129,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4134712\nallocs=10990\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 24111712.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=29248\nallocs=394\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "f5b9d1a3516160a536e0be29a0adc447743c2e95",
          "message": "Fix ring_collision test",
          "timestamp": "2026-03-19T15:21:29+01:00",
          "tree_id": "47649aa3c54f5e52a3f664c758c6f427bbd1b3be",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/f5b9d1a3516160a536e0be29a0adc447743c2e95"
        },
        "date": 1773930698059,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/KaiserBesselSplitting/add_local_integrals!",
            "value": 12221796.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22368\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting/velocity + streamfunction",
            "value": 562382417,
            "unit": "ns",
            "extra": "gctime=0\nmemory=230544\nallocs=2294\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting/velocity",
            "value": 441532221,
            "unit": "ns",
            "extra": "gctime=0\nmemory=183488\nallocs=1801\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 12173244,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20320\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12464440.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5088\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 704983480.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=226848\nallocs=2294\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 586638060,
            "unit": "ns",
            "extra": "gctime=0\nmemory=179968\nallocs=1801\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/add_local_integrals!",
            "value": 12181954,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20320\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/velocity + streamfunction",
            "value": 707514747,
            "unit": "ns",
            "extra": "gctime=0\nmemory=226848\nallocs=2294\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/velocity",
            "value": 591484077,
            "unit": "ns",
            "extra": "gctime=0\nmemory=179968\nallocs=1801\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8700905,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1150301916,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4690296\nallocs=7516\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1607507836.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42130256\nallocs=25352\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5461576,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 813346,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2236316,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 4098982708,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68636176\nallocs=176606\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 219787276,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9931848\nallocs=2697\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 7068564,
            "unit": "ns",
            "extra": "gctime=0\nmemory=99376\nallocs=2593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 2120344.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (sorted)",
            "value": 55267497,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/unsorted)",
            "value": 85402399,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 75166279.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 71113178.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/sorted)",
            "value": 78333710,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 61596848,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 5777454,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (sorted)",
            "value": 93601549.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/unsorted)",
            "value": 256689381,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 338845658,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 258729356.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/sorted)",
            "value": 107574558,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 260920453.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2740952594,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4134712\nallocs=10990\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 24674087,
            "unit": "ns",
            "extra": "gctime=0\nmemory=29248\nallocs=394\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "49699333+dependabot[bot]@users.noreply.github.com",
            "name": "dependabot[bot]",
            "username": "dependabot[bot]"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "a1fa03594059be7b3f3c73179f6f9c060e139c25",
          "message": "Bump codecov/codecov-action from 5 to 6 (#98)\n\nBumps [codecov/codecov-action](https://github.com/codecov/codecov-action) from 5 to 6.\n- [Release notes](https://github.com/codecov/codecov-action/releases)\n- [Changelog](https://github.com/codecov/codecov-action/blob/main/CHANGELOG.md)\n- [Commits](https://github.com/codecov/codecov-action/compare/v5...v6)\n\n---\nupdated-dependencies:\n- dependency-name: codecov/codecov-action\n  dependency-version: '6'\n  dependency-type: direct:production\n  update-type: version-update:semver-major\n...\n\nSigned-off-by: dependabot[bot] <support@github.com>\nCo-authored-by: dependabot[bot] <49699333+dependabot[bot]@users.noreply.github.com>",
          "timestamp": "2026-03-30T11:34:43+02:00",
          "tree_id": "ccdec8bb8c08ee3b294f1728f0a50432cf4a5c4c",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/a1fa03594059be7b3f3c73179f6f9c060e139c25"
        },
        "date": 1774864417731,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/KaiserBesselSplitting/add_local_integrals!",
            "value": 12535346,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22368\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting/velocity + streamfunction",
            "value": 585683871,
            "unit": "ns",
            "extra": "gctime=0\nmemory=236592\nallocs=2396\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting/velocity",
            "value": 474490652,
            "unit": "ns",
            "extra": "gctime=0\nmemory=189536\nallocs=1903\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 12508261,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20320\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12873968,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5088\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 743304623,
            "unit": "ns",
            "extra": "gctime=0\nmemory=232896\nallocs=2396\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 619156190,
            "unit": "ns",
            "extra": "gctime=0\nmemory=186016\nallocs=1903\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/add_local_integrals!",
            "value": 12506697.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20320\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/velocity + streamfunction",
            "value": 745496621,
            "unit": "ns",
            "extra": "gctime=0\nmemory=232896\nallocs=2396\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/velocity",
            "value": 625371265.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=186016\nallocs=1903\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8839305,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1212678335,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4690296\nallocs=7516\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1601778454.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42130256\nallocs=25352\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5569669,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 820650,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2248999,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 4001303257.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68636432\nallocs=176606\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 220959135,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9931848\nallocs=2697\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 7079654,
            "unit": "ns",
            "extra": "gctime=0\nmemory=99376\nallocs=2593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 2113962,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (sorted)",
            "value": 64446106,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/unsorted)",
            "value": 115811445,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 101900712,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 83755819,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/sorted)",
            "value": 87388534,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 82243948.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 7699349,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (sorted)",
            "value": 111550115,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/unsorted)",
            "value": 307951688,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 396248070,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 300919763,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/sorted)",
            "value": 123666669,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 308109474,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2857155542,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4158792\nallocs=11397\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 23896286,
            "unit": "ns",
            "extra": "gctime=0\nmemory=29248\nallocs=394\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "3abdb6833daca8695c81e00e29d1af5a40f65714",
          "message": "Update rotation example",
          "timestamp": "2026-04-01T15:58:51+02:00",
          "tree_id": "cabc026bc134987032a08ddf7fd99714c6788b4d",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/3abdb6833daca8695c81e00e29d1af5a40f65714"
        },
        "date": 1775053010442,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/KaiserBesselSplitting/add_local_integrals!",
            "value": 12221915,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22368\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting/velocity + streamfunction",
            "value": 570887408,
            "unit": "ns",
            "extra": "gctime=0\nmemory=230544\nallocs=2294\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting/velocity",
            "value": 450896307,
            "unit": "ns",
            "extra": "gctime=0\nmemory=183488\nallocs=1801\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 12243545,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20320\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12577959,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5088\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 720970916,
            "unit": "ns",
            "extra": "gctime=0\nmemory=226848\nallocs=2294\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 598977779,
            "unit": "ns",
            "extra": "gctime=0\nmemory=179968\nallocs=1801\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/add_local_integrals!",
            "value": 12292343,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20320\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/velocity + streamfunction",
            "value": 727118972,
            "unit": "ns",
            "extra": "gctime=0\nmemory=226848\nallocs=2294\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/velocity",
            "value": 608302368,
            "unit": "ns",
            "extra": "gctime=0\nmemory=179968\nallocs=1801\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8681632.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1193968329,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4690296\nallocs=7516\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1535797090.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42130256\nallocs=25352\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5497639,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 817509,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2238405,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3921438559.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68636176\nallocs=176606\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 224367715,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9931848\nallocs=2697\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 6974301,
            "unit": "ns",
            "extra": "gctime=0\nmemory=99376\nallocs=2593\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 2104847,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (sorted)",
            "value": 56989722,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/unsorted)",
            "value": 84670636,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 75789964,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 62168854.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/sorted)",
            "value": 76862968,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 65998424.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 6538826,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (sorted)",
            "value": 103984146,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/unsorted)",
            "value": 283434475.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 353311291,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 267437815,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/sorted)",
            "value": 115718611,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 274366582,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2806943715,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4134712\nallocs=10990\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10,\"evals\":1,\"gcsample\":false,\"seconds\":1000,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 23989510.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=29248\nallocs=394\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "20c82e952fe26f18f35753e657e1736824eb673e",
          "message": "v0.33.2",
          "timestamp": "2026-04-13T13:20:07+02:00",
          "tree_id": "eecbd455c6590ab36c8da44e4c7bd63da5b73097",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/20c82e952fe26f18f35753e657e1736824eb673e"
        },
        "date": 1776080023293,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/GaussianSplitting/add_local_integrals!",
            "value": 13464462,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20320\nallocs=34\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/velocity",
            "value": 654297164,
            "unit": "ns",
            "extra": "gctime=0\nmemory=186016\nallocs=1903\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/velocity + streamfunction",
            "value": 776550938,
            "unit": "ns",
            "extra": "gctime=0\nmemory=238944\nallocs=2498\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting/add_local_integrals!",
            "value": 13473206,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22368\nallocs=34\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting/velocity",
            "value": 487148470,
            "unit": "ns",
            "extra": "gctime=0\nmemory=189536\nallocs=1903\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting/velocity + streamfunction",
            "value": 608239950,
            "unit": "ns",
            "extra": "gctime=0\nmemory=242640\nallocs=2498\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 13467280,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20320\nallocs=34\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 13662874,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5088\nallocs=34\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 652520228.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=186016\nallocs=1903\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 775969933,
            "unit": "ns",
            "extra": "gctime=0\nmemory=238944\nallocs=2498\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/sorted)",
            "value": 85591317,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/unsorted)",
            "value": 92363771,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (sorted)",
            "value": 61374839.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 68192710,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 67927505,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 84149740,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 2101618,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/sorted)",
            "value": 116770501,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/unsorted)",
            "value": 230838971.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (sorted)",
            "value": 114246440.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 300966865,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 310488108,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 381421786,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 6401432,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1247329716,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4744728\nallocs=8434\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 9602862.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2353833,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1608142558,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42184688\nallocs=26270\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 6152657,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 916249,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 4398955088,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68636112\nallocs=176606\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 188459091.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9931848\nallocs=2697\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 7043575,
            "unit": "ns",
            "extra": "gctime=0\nmemory=99376\nallocs=2593\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 25657781,
            "unit": "ns",
            "extra": "gctime=0\nmemory=29248\nallocs=394\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Timestepping/step!",
            "value": 3035702574.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4164952\nallocs=11500\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10,\"seconds\":1000,\"time_tolerance\":0.05}"
          }
        ]
      },
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
          "id": "f28ebc1dcfb993e6b88b888d307f818bca4de07b",
          "message": "Disable OpenCLBackend tests on CI",
          "timestamp": "2026-04-14T10:13:17+02:00",
          "tree_id": "7bff25842b64523a4d175932ec380e5702b6b702",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/f28ebc1dcfb993e6b88b888d307f818bca4de07b"
        },
        "date": 1776154989961,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/GaussianSplitting/add_local_integrals!",
            "value": 12453556,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20320\nallocs=34\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/velocity",
            "value": 628639276.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=192064\nallocs=2005\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/velocity + streamfunction",
            "value": 757054767,
            "unit": "ns",
            "extra": "gctime=0\nmemory=244992\nallocs=2600\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting/add_local_integrals!",
            "value": 12460852.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22368\nallocs=34\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting/velocity",
            "value": 471446934,
            "unit": "ns",
            "extra": "gctime=0\nmemory=195584\nallocs=2005\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting/velocity + streamfunction",
            "value": 597069403,
            "unit": "ns",
            "extra": "gctime=0\nmemory=248688\nallocs=2600\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 12416200.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20320\nallocs=34\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 12766193,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5088\nallocs=34\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 610883019,
            "unit": "ns",
            "extra": "gctime=0\nmemory=192064\nallocs=2005\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 737621077,
            "unit": "ns",
            "extra": "gctime=0\nmemory=244992\nallocs=2600\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/sorted)",
            "value": 88548129,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/unsorted)",
            "value": 107210586.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (sorted)",
            "value": 61964992,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 74258813,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 74504653,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 90917867,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 2102992.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/sorted)",
            "value": 128201847.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/unsorted)",
            "value": 303292967,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (sorted)",
            "value": 105372184,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 302409144,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 305110108,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 394441964,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 8093034,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1225220014,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4744728\nallocs=8434\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 8793045,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2265929,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1521708625,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42184688\nallocs=26270\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 5511705,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 817774,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 3984586484.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68636432\nallocs=176606\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 236651149,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9931848\nallocs=2697\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 6924579,
            "unit": "ns",
            "extra": "gctime=0\nmemory=99376\nallocs=2593\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 24668839,
            "unit": "ns",
            "extra": "gctime=0\nmemory=29248\nallocs=394\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2864216186,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4189032\nallocs=11907\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10,\"seconds\":1000,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "49699333+dependabot[bot]@users.noreply.github.com",
            "name": "dependabot[bot]",
            "username": "dependabot[bot]"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "2c331babb54d8b5f8d05e6432591a4eae68bc5a9",
          "message": "Bump julia-actions/setup-julia from 2 to 3 (#100)\n\nBumps [julia-actions/setup-julia](https://github.com/julia-actions/setup-julia) from 2 to 3.\n- [Release notes](https://github.com/julia-actions/setup-julia/releases)\n- [Commits](https://github.com/julia-actions/setup-julia/compare/v2...v3)\n\n---\nupdated-dependencies:\n- dependency-name: julia-actions/setup-julia\n  dependency-version: '3'\n  dependency-type: direct:production\n  update-type: version-update:semver-major\n...\n\nSigned-off-by: dependabot[bot] <support@github.com>\nCo-authored-by: dependabot[bot] <49699333+dependabot[bot]@users.noreply.github.com>",
          "timestamp": "2026-04-20T10:51:55+02:00",
          "tree_id": "9eba47c7d806c0b46d3bf6e28311629a51ce450c",
          "url": "https://github.com/jipolanco/VortexPasta.jl/commit/2c331babb54d8b5f8d05e6432591a4eae68bc5a9"
        },
        "date": 1776675720400,
        "tool": "julia",
        "benches": [
          {
            "name": "BiotSavart/GaussianSplitting/add_local_integrals!",
            "value": 13564538,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20320\nallocs=34\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/velocity",
            "value": 670304551,
            "unit": "ns",
            "extra": "gctime=0\nmemory=186016\nallocs=1903\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "BiotSavart/GaussianSplitting/velocity + streamfunction",
            "value": 793395834,
            "unit": "ns",
            "extra": "gctime=0\nmemory=238944\nallocs=2498\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting/add_local_integrals!",
            "value": 13558860,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22368\nallocs=34\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting/velocity",
            "value": 494912756,
            "unit": "ns",
            "extra": "gctime=0\nmemory=189536\nallocs=1903\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "BiotSavart/KaiserBesselSplitting/velocity + streamfunction",
            "value": 622440772,
            "unit": "ns",
            "extra": "gctime=0\nmemory=242640\nallocs=2498\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "BiotSavart/add_local_integrals!",
            "value": 13475495,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20320\nallocs=34\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "BiotSavart/add_point_charges!",
            "value": 13845374,
            "unit": "ns",
            "extra": "gctime=0\nmemory=5088\nallocs=34\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "BiotSavart/velocity",
            "value": 663659895.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=186016\nallocs=1903\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "BiotSavart/velocity + streamfunction",
            "value": 786769657,
            "unit": "ns",
            "extra": "gctime=0\nmemory=238944\nallocs=2498\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/sorted)",
            "value": 87239152,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (SIMD/unsorted)",
            "value": 97823376,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (sorted)",
            "value": 64587045.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_pair (unsorted)",
            "value": 70805416.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/foreach_source",
            "value": 72360670,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/iterator_interface",
            "value": 84786537,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 1/set_elements!",
            "value": 2133463.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/sorted)",
            "value": 117014234,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (SIMD/unsorted)",
            "value": 234201151,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (sorted)",
            "value": 120218356,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6112\nallocs=66\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_pair (unsorted)",
            "value": 295244609,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=22\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/foreach_source",
            "value": 280458049,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/iterator_interface",
            "value": 352307396,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2208\nallocs=22\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "CellLists/CPU/nsubdiv = 2/set_elements!",
            "value": 6571460,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3776\nallocs=44\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Diagnostics/energy_flux",
            "value": 1217263114,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4744728\nallocs=8434\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Diagnostics/energy_injection_rate",
            "value": 9616654,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Diagnostics/energy_spectrum",
            "value": 2373617.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3296\nallocs=26\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Diagnostics/energy_transfer_matrix",
            "value": 1579241600,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42184688\nallocs=26270\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Diagnostics/helicity",
            "value": 6158036,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Diagnostics/kinetic_energy",
            "value": 917112,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4528\nallocs=39\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Reconnections/ReconnectBasedOnDistance",
            "value": 4328081316,
            "unit": "ns",
            "extra": "gctime=0\nmemory=68636432\nallocs=176606\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Reconnections/ReconnectFast",
            "value": 190937364.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9931848\nallocs=2697\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Refinement/RefineBasedOnSegmentLength",
            "value": 7225481.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=99376\nallocs=2593\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Timestepping/forcing",
            "value": 26542652,
            "unit": "ns",
            "extra": "gctime=0\nmemory=29248\nallocs=394\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Timestepping/step!",
            "value": 2929188379.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=4164952\nallocs=11500\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10,\"seconds\":1000,\"time_tolerance\":0.05}"
          }
        ]
      }
    ]
  }
}