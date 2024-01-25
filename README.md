# ticinfo

I wanted a way to quickly see some info about a star based on it's TICID, so I made this tool.

Maybe someone else will find it useful

```
pip install ticinfo
```

to use

```bash
% toco 441420236

    ID        ra       dec     pmRA  pmDEC    eclong     eclat    Tmag Vmag Kmag Teff   rad      mass     d
--------- ---------- -------- ------ ------ ---------- ---------- ---- ---- ---- ---- -------- -------- ------
441420236 311.289718 -31.3409 281.42 -359.9 305.309919 -12.822828 6.76 8.81 4.53  nan 0.698009 0.662074 9.7221

Target name: V* AU Mic
The target is in constellation Microscopium
```

if you just have coordiantes you can use

```bash
% toco 311.289718 -31.3409
    ID        ra       dec     pmRA  pmDEC    eclong     eclat    Tmag Vmag Kmag Teff   rad      mass     d
--------- ---------- -------- ------ ------ ---------- ---------- ---- ---- ---- ---- -------- -------- ------
441420236 311.289718 -31.3409 281.42 -359.9 305.309919 -12.822828 6.76 8.81 4.53  nan 0.698009 0.662074 9.7221

Target name: V* AU Mic
The target is in constellation Microscopium
```

or, if you know the name of the source you can use

```bash
% toco au mic
     ID        ra       dec     pmRA  pmDEC    eclong     eclat    Tmag Vmag Kmag Teff   rad      mass     d
--------- ---------- -------- ------ ------ ---------- ---------- ---- ---- ---- ---- -------- -------- ------
441420236 311.289718 -31.3409 281.42 -359.9 305.309919 -12.822828 6.76 8.81 4.53  nan 0.698009 0.662074 9.7221

Target name: V* AU Mic
The target is in constellation Microscopium
```

## Finding whether data is in the archive

You are able to find out whether there is TESS data available in the archive for any target by using the `--stored` flag, or `-s` for short.

```bash
% toco L 98-59 -s
    ID        ra       dec    pmRA  pmDEC    eclong     eclat    Tmag  Vmag Kmag  Teff    rad     mass      d
--------- ---------- ------- ----- ------- ---------- ---------- ---- ----- ---- ------ ------- -------- -------
307210830 124.531756 -68.313 94.77 -340.47 203.322782 -76.813828 9.41 11.68  7.1 3429.0 0.31416 0.292836 10.6194

Target name: L   98-59
The target is in constellation Volans
FFI data at MAST for sectors:   [1, 2, 5, 8, 9, 10, 11, 12]
2-min data at MAST for sectors: [2, 5, 8, 9, 10, 11, 12]
20-s data at MAST for sectors:  []
```

## Finding whether data will be available

You are able to find out whether there the target will be observable in the next TESS cycle using the `--future` flag, or `-f` for short.

```bash
% toco L 98-59 -f
The target is in constellation Volans                                               
                                                                                    
    ID        ra      dec    plx  Tmag  Vmag Kmag  Teff    rad     mass             
--------- --------- ------- ----- ---- ----- ---- ------ ------- --------           
307210830 124.53176 -68.313 94.14 9.41 11.68  7.1 3429.0 0.31416 0.292836           
                                                                                    
Future: Will be observable in Cycle 7 sectors:  [88 89 90 96]                       
```
