
Simple example of Kirchhof seismic migration.

## usage:

```
make cython
make run
```

#### Original data

![](./figures/unprocessed.png)

#### Migrated data

![](./figures/migrated_data.png)

There are some artifacts visible, which we can improve by tapering the edges:

#### Tapered original data

![](./figures/tapered_original_data.png)

#### Migrated taperad data

![](./figures/migrated_tapered_data.png)
