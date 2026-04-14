# Updating

## Updating databases

### Rerun PointFinder database creation
```
rm assets/pointfinder_db/*/*.{b,name}
make update_pointfinder_db
```

### Rerun SerotypeFinder database creation
```
rm assets/serotypefinder_db/*.{b,name}
make update_serotypefinder_db
```

### Rerun ResFinder database creation
```
rm assets/resfinder_db/*.{b,name}
make update_resfinder_db
```

### Rerun VirulenceFinder database creation
```
rm assets/virulencefinder_db/*.{b,name}
make update_virulencefinder_db
```

### Update AMRFinder database
```
rm assets/amrfinder_db/*.{b,name}
make update_amrfinderplus
```