# Updating

## Submodules

### Simultaneously update all submodules to expected commit
```
git submodule update --init --recursive
```

### Check that everything worked
**NOTE**: Submodules should not appear when running the following
```
git status
```

### Troubleshooting
Git submodule update not working
```
cd assets/{submodule_dir}
git pull origin master
git checkout commit_id
cd ../..
git status
```

### Rerun PointFinder database creation
```
rm assets/pointfinder_db/*/*.{b,name}
make
```

### Rerun SerotypeFinder database creation
```
rm assets/serotypefinder_db/*.{b,name}
make
```

### Rerun ResFinder database creation
```
rm assets/resfinder_db/*.{b,name}
make
```

### Rerun VirulenceFinder database creation
```
rm assets/virulencefinder_db/*.{b,name}
make
```
