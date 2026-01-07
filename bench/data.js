window.BENCHMARK_DATA = {
  "lastUpdate": 1767780870766,
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
      }
    ]
  }
}