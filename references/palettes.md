# IOBR Palette and Theme Reference

## design_mytheme()

IOBR provides a unified theme function `design_mytheme()` for consistent Nature-style publication figures.
Use this for all ggplot-based visualizations.

## Categorical Palettes (`palette_group` parameter)

Use these for grouping variables (tissue types, TME subtypes, clinical groups).

### Journal Palettes
| Name | Colors | Best For |
|------|--------|----------|
| `"jama"` | 7 colors, bold primary | Default choice, clinical studies |
| `"aaas"` | 10 colors, Science journal | General publication |
| `"jco"` | 10 colors, clinical oncology | Oncology papers |
| `"nrc"` | Nature Research colors | Nature journals |

### Functional Palettes
| Name | Colors | Best For |
|------|--------|----------|
| `"paired1"` | 6 paired colors | 2-condition comparisons |
| `"paired2"` | 6 paired colors | Alternative pairing |
| `"paired3"` | 6 paired colors | Alternative pairing |
| `"paired4"` | 6 paired colors | Alternative pairing |
| `"accent"` | 8 colors | Colorblind-friendly |
| `"set2"` | 8 colors, soft | Many groups, gentle appearance |

## Heatmap Palettes (`palette` parameter, numeric 1–6)

Use these for continuous data in heatmaps (expression, signature scores).

| Number | Color Range | Description |
|--------|------------|-------------|
| 1 | Classic pheatmap | Red-blue traditional |
| 2 | Blue → White → Red | Diverging, default |
| 3 | Purple → White → Orange | Warm diverging |
| 4 | White → Blue | Sequential blue |
| 5 | White → Red/Orange | Sequential warm |
| 6 | White → Teal/Cyan | Sequential cool |

## Custom Colors

Both `sig_heatmap()` and other plot functions accept custom colors:

```r
# Custom group colors (comma-separated hex)
sig_heatmap(..., custom_cols = "#E64B35,#4DBBD5,#00A087")

# Custom heatmap gradient (at least 3 colors)
sig_heatmap(..., custom_heatmap_cols = c("#2166AC", "#F7F7F7", "#B2182B"))
```

## Selection Guide

When unsure, use this decision tree:

1. **How many groups?**
   - 2–3 groups → `"jama"` or `"jco"`
   - 4–6 groups → `"aaas"` or `"nrc"`
   - 7+ groups → `"set2"` or `"accent"`
2. **Is it a heatmap?** → Start with palette `2` (blue-white-red)
3. **Need colorblind safety?** → `"accent"` or `"set2"`
4. **Clinical/oncology context?** → `"jco"` or `"jama"`
5. **Nature publication?** → `"nrc"`

If still unsure, present 3 options to the user and let them choose.
