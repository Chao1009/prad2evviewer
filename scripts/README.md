# Scripts

Utility scripts for PRad-II detector visualization and analysis. Requires Python 3 with matplotlib and numpy.

```bash
pip install -r scripts/requirements.txt
```

## gem_layout.py

Visualize GEM strip layout from `gem_map.json`. Shows X/Y strips and APV boundaries with beam hole.

```bash
python scripts/gem_layout.py [path/to/gem_map.json]
```

Defaults to `database/gem_map.json`. Saves `gem_layout.png`.

## gem_cluster_view.py

Visualize GEM clustering results from a `gem_dump -m evdump` JSON file. Overlays fired strips, clusters, and 2D hits on the detector geometry.

```bash
python scripts/gem_cluster_view.py <event.json> [gem_map.json] [--det N] [-o file.png]
```

| Option | Description |
|--------|-------------|
| `event.json` | Event dump from `gem_dump -m evdump` (required) |
| `gem_map.json` | GEM map file (default: auto-search `database/gem_map.json`) |
| `--det N` | Show only detector N (default: all detectors) |
| `-o file` | Output image (default: `gem_cluster_view.png`) |

Visual elements:
- **X strip hits** — vertical lines, blue colormap scaled by charge
- **Y strip hits** — horizontal lines, red colormap scaled by charge
- **Cross-talk hits** — dashed lines at lower opacity
- **Cluster ranges** — semi-transparent bands over constituent strips
- **Cluster centers** — triangle markers at detector edge with charge/size labels
- **2D hits** — green star markers at reconstructed (x, y) positions
- **Beam hole** — yellow rectangle

Example workflow:
```bash
# 1. Compute pedestals (if not already done)
gem_dump data.evio -m ped -o gem_ped.json

# 2. Dump a single event to JSON
gem_dump data.evio -P gem_ped.json -m evdump -e 42 -o event_42.json

# 3. Visualize clustering
python scripts/gem_cluster_view.py event_42.json database/gem_map.json

# Show only GEM0
python scripts/gem_cluster_view.py event_42.json --det 0 -o gem0_event42.png
```
