# UTIA BRDF Library
Python script to read data from the <a href="http://btf.utia.cas.cz/?brdf_dat_dwn">UTIA BRDF database</a><sup>1</sup>. The database consists of 150 files in binary, MATLAB, and OpenEXR formats (color-space sRGB, observer 2 degrees, illuminant D65).

<br>

<div align="center">

![spheres](https://user-images.githubusercontent.com/10238412/127757059-9475f99f-b512-4452-ad00-12be46f75eac.jpg)

</div>

### Usage
Download BRDF data from the <a href="http://btf.utia.cas.cz/?brdf_dat_dwn">UTIA BRDF database</a> in .bin format into the _data_ folder and run:
```$ python read_brdf.py```

### References
[1] Template-Based Sampling of Anisotropic BRDFs (Filip J., Vavra R.),  Computer Graphics Forum (Proceedings of Pacific Graphics 2014, Seoul, Korea), Eurographics 2014.
