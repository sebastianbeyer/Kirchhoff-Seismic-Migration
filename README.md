
Simple implementation of [Kirchhof seismic migration](https://en.wikipedia.org/wiki/Seismic_migration) in Python/Cython.

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
