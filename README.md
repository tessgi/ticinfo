# toco
I wanted a way to quickly see some info about a star based on it's TICID, so I made this tool.

Maybe someone else will find it useful

```
pip install ticinfo
```

to use

```bash
% toco 471011144

    ID        ra        dec       pmRA   pmDEC    eclong     eclat    Tmag Vmag Kmag Teff rad mass    d
--------- ---------- ---------- -------- ------ ---------- --------- ----- ---- ---- ---- --- ---- --------
471011144 219.903981 -60.837157 -3600.35 952.11 239.482264 -42.59683 -0.15 1.35   --   --  --   -- 1.347491
```


if you have a star you think might be in Simbad, you can do this
```bash
% tocot 471011144

    ID        ra        dec       pmRA   pmDEC    eclong     eclat    Tmag Vmag Kmag Teff rad mass    d
--------- ---------- ---------- -------- ------ ---------- --------- ----- ---- ---- ---- --- ---- --------
471011144 219.903981 -60.837157 -3600.35 952.11 239.482264 -42.59683 -0.15 1.35   --   --  --   -- 1.347491

Target name: * alf Cen
The target is in constellation Centaurus
```

if you jusat have coordiantes you can use 
```bash
% tococ 219.90085 -60.83561944444
    ID        ra        dec       pmRA   pmDEC    eclong     eclat    Tmag Vmag Kmag Teff rad mass    d
--------- ---------- ---------- -------- ------ ---------- --------- ----- ---- ---- ---- --- ---- --------
471011144 219.903981 -60.837157 -3600.35 952.11 239.482264 -42.59683 -0.15 1.35   --   --  --   -- 1.347491

Target name: * alf Cen
The target is in constellation Centaurus
```

or, if you know the name of the source you can use
```bash
% tocon alpha cen
    ID        ra        dec       pmRA   pmDEC    eclong     eclat    Tmag Vmag Kmag Teff rad mass    d
--------- ---------- ---------- -------- ------ ---------- --------- ----- ---- ---- ---- --- ---- --------
471011144 219.903981 -60.837157 -3600.35 952.11 239.482264 -42.59683 -0.15 1.35   --   --  --   -- 1.347491

Target name: * alf Cen
The target is in constellation Centaurus
```