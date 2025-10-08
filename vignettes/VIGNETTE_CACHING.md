# Vignette Caching System

This package uses an intelligent vignette caching system to dramatically reduce build times while maintaining the ability to force full rebuilds when needed.

## How It Works

### Default Behavior (Fast Builds)
- **Vignettes default to using cached data** (`runchunks = FALSE`)
- **Build times are ~90% faster** when using cached data
- **Each vignette has its own cache file** for independent caching

### Force Rebuild Options

#### 1. Environment Variable
```bash
# Force rebuild all vignettes
FORCE_VIGNETTE_REBUILD=TRUE Rscript -e "tools::buildVignettes(dir = '.', tangle = TRUE)"

# Use cached data (default)
Rscript -e "tools::buildVignettes(dir = '.', tangle = TRUE)"
```

#### 2. Build Script
```bash
# Fast build (use cache)
./build_vignettes.R

# Force rebuild
./build_vignettes.R --force
```

#### 3. GitHub Actions
```bash
# Force rebuild in CI (add to commit message)
git commit -m "Update vignettes /rebuild-vignettes"
```

## Cache Files

Each vignette has its own cache file:
- `vignettes/SpiecEasi.RData` - Main vignette
- `vignettes/phyloseq-integration.RData` - Phyloseq integration
- `vignettes/cross-domain-interactions.RData` - Cross-domain analysis
- `vignettes/latent-variable-models.RData` - Latent variable models
- `vignettes/pulsar-parallel.RData` - Parallel processing
- `vignettes/troubleshooting.RData` - Troubleshooting

## When to Force Rebuild

Force a full rebuild when:
- ✅ **Vignette content changes** (R code, text, parameters)
- ✅ **Package functions change** that affect vignette results
- ✅ **Data processing changes** that affect cached objects
- ✅ **First-time builds** (no cache exists)

Use cached data when:
- ✅ **Only documentation changes** (text, formatting)
- ✅ **No functional changes** to vignette code
- ✅ **Quick builds for testing**

## Performance Benefits

| Build Type | Time | Use Case |
|------------|------|----------|
| **Cached** | ~2-3 minutes | Documentation updates, quick testing |
| **Full Rebuild** | ~15-20 minutes | Code changes, first build |

## Troubleshooting

### Cache Issues
```bash
# Clear all caches and force rebuild
rm vignettes/*.RData
./build_vignettes.R --force
```

### Individual Vignette
```bash
# Rebuild specific vignette
rm vignettes/SpiecEasi.RData
FORCE_VIGNETTE_REBUILD=TRUE Rscript -e "rmarkdown::render('vignettes/SpiecEasi.Rmd')"
```

### CI/CD Integration
The GitHub Actions workflow automatically:
- ✅ **Caches individual vignette files**
- ✅ **Detects `/rebuild-vignettes` in commit messages**
- ✅ **Falls back to cached builds by default**
- ✅ **Provides detailed logging of cache usage**
