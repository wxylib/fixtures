"""Index space against coordinate space (design m07, section 5): small
Zarr v3 stores as a real producer writes them, plus the truth that
producer threw away, so a reader can be checked on both.

I1 roms-fan: a ROMS-style coastal grid, regular in index space and
curvilinear on the globe (playproj's `fan()`, log-polar in stereographic
radius), with a terrain-following s-coordinate. Stored with ROMS's own
names and CF attributes: `lon_rho`, `lat_rho`, `h`, `zeta`, `mask_rho`,
`s_rho` (CF `ocean_s_coordinate_g2` with `formula_terms`), `Cs_r`, `hc`,
`Vtransform`, `Vstretching`, `theta_s`, `theta_b` and
`temp(s_rho, eta_rho, xi_rho)`. Oracles: depths of every rho and w level
from a NumPy port of ROMS set_scoord.F and set_depth.F (the one in
playproj's test/make-roms-fixture.py), and GDAL's geolocation-array warp
of the top layer to EPSG:4326 (osgeo.gdal).

I3 rotated: a coarse EURO-CORDEX window on a CF
`rotated_latitude_longitude` grid: 1D `rlon`, `rlat`, the grid-mapping
variable `rotated_pole`, and float32 2D `lon`, `lat`, as CORDEX files
carry them. The grid is regular (R1) in PROJ's `ob_tran`; the manifest
gives the generating affine and the PROJ string, and osgeo.osr checks
the stored 2D arrays against them.

Run from the repository root, with a Python that has GDAL's bindings,
numpy, xarray and zarr (3.x):

    python scripts/make_index_fixtures.py
"""

import json
import os
import shutil

import numpy as np
import osgeo
import xarray as xr
from osgeo import gdal, osr

gdal.UseExceptions()

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "index")
D2R = np.pi / 180.0
METHODS = ("near", "bilinear")


def wrap(lon):
    return (lon + 180.0) % 360.0 - 180.0


def destination(lon0, lat0, c, az):
    """Point at angular distance c (radians) and azimuth az from lon0, lat0
    on the sphere."""
    p0, l0 = lat0 * D2R, lon0 * D2R
    lat = np.arcsin(np.sin(p0) * np.cos(c) + np.cos(p0) * np.sin(c) * np.cos(az))
    lon = l0 + np.arctan2(np.sin(az) * np.sin(c) * np.cos(p0),
                          np.cos(c) - np.sin(p0) * np.sin(lat))
    return wrap(lon / D2R), lat / D2R


# ---------------------------------------------------------------- I1

def fan(ni=24, lon0=150.4, lat0=-33.2, bearing=105, spread=90, r0=0.9, r1=18):
    """playproj js/roms.js fan(), at cell centres: equal steps in bearing
    across the fan (xi) and in log stereographic radius along it (eta).
    Returns rho-point lon, lat shaped (eta, xi)."""
    q0 = np.log(2 * np.tan(r0 * D2R / 2))
    q1 = np.log(2 * np.tan(r1 * D2R / 2))
    a = spread * D2R
    nj = max(4, round(ni * (q1 - q0) / a))
    u = (np.arange(ni) + 0.5) / ni
    v = (np.arange(nj) + 0.5) / nj
    U, V = np.meshgrid(u, v)
    c = 2 * np.arctan(np.exp(q0 + (q1 - q0) * V) / 2)
    return destination(lon0, lat0, c, bearing * D2R + (U - 0.5) * a)


# ROMS set_scoord.F and set_depth.F, ported to NumPy in playproj's
# test/make-roms-fixture.py (copied here unchanged in substance).

def cs(s, Vstretching, theta_s, theta_b):
    s = np.asarray(s, dtype=float)
    if Vstretching == 1:
        if theta_s > 0:
            ptheta = np.sinh(theta_s * s) / np.sinh(theta_s)
            rtheta = np.tanh(theta_s * (s + 0.5)) / (2.0 * np.tanh(0.5 * theta_s)) - 0.5
            return (1.0 - theta_b) * ptheta + theta_b * rtheta
        return s
    if Vstretching == 2:
        Aweight, Bweight = 1.0, 1.0
        if theta_s > 0:
            Csur = (1.0 - np.cosh(theta_s * s)) / (np.cosh(theta_s) - 1.0)
            if theta_b > 0:
                Cbot = -1.0 + np.sinh(theta_b * (s + 1.0)) / np.sinh(theta_b)
                Cweight = (s + 1.0) ** Aweight * (1.0 + (Aweight / Bweight) * (1.0 - (s + 1.0) ** Bweight))
                return Cweight * Csur + (1.0 - Cweight) * Cbot
            return Csur
        return s
    if Vstretching == 4:
        if theta_s > 0:
            Csur = (1.0 - np.cosh(theta_s * s)) / (np.cosh(theta_s) - 1.0)
        else:
            Csur = -s ** 2
        if theta_b > 0:
            return (np.exp(theta_b * Csur) - 1.0) / (1.0 - np.exp(-theta_b))
        return Csur
    raise ValueError(Vstretching)


def depth(s, C, h, zeta, Vtransform, hc):
    """z at levels s (stretching C) over h and zeta: s on the first axis."""
    s = np.asarray(s)[:, None, None]
    C = np.asarray(C)[:, None, None]
    if Vtransform == 1:
        z0 = hc * s + (h - hc) * C
        return z0 + zeta * (1.0 + z0 / h)
    z0 = (hc * s + h * C) / (hc + h)
    return zeta + (zeta + h) * z0


VERT = {"Vtransform": 2, "Vstretching": 4, "theta_s": 6.5, "theta_b": 1.5, "hc": 200.0}
N = 8


def write_tif(path, bands, dtype):
    ny, nx = bands[0].shape
    d = gdal.GetDriverByName("GTiff").Create(path, nx, ny, len(bands), dtype)
    for k, a in enumerate(bands):
        d.GetRasterBand(k + 1).WriteArray(a)
    d = None


def write_geoloc_vrt(path, src_tif, lonlat_tif, ny, nx):
    srs = osr.SpatialReference()
    srs.SetFromUserInput("EPSG:4326")
    md = {
        "X_DATASET": lonlat_tif, "X_BAND": "1",
        "Y_DATASET": lonlat_tif, "Y_BAND": "2",
        "PIXEL_OFFSET": "0", "LINE_OFFSET": "0",
        "PIXEL_STEP": "1", "LINE_STEP": "1",
        "SRS": srs.ExportToWkt(),
        "GEOREFERENCING_CONVENTION": "PIXEL_CENTER",
    }
    esc = lambda v: v.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")
    items = "\n".join('    <MDI key="%s">%s</MDI>' % (k, esc(v)) for k, v in md.items())
    with open(path, "w") as f:
        f.write("""<VRTDataset rasterXSize="%d" rasterYSize="%d">
  <Metadata domain="GEOLOCATION">
%s
  </Metadata>
  <VRTRasterBand dataType="Float32" band="1">
    <SimpleSource>
      <SourceFilename relativeToVRT="1">%s</SourceFilename>
      <SourceBand>1</SourceBand>
    </SimpleSource>
  </VRTRasterBand>
</VRTDataset>
""" % (nx, ny, items, os.path.basename(src_tif)))


def roms_fan():
    lon, lat = fan()
    ny, nx = lon.shape
    jj, ii = np.mgrid[0:ny, 0:nx]
    # shelf near the coast (eta 0), deepening offshore, a canyon across xi
    h = 15.0 + 3000.0 / (1.0 + np.exp(-(jj - 0.45 * ny) / 4.0)) \
        + 400.0 * np.exp(-((ii - 0.6 * nx) / 2.5) ** 2) * (jj > 4)
    zeta = 0.3 * np.sin(ii / nx * np.pi) * np.exp(-jj / 10.0)
    mask = np.ones((ny, nx))
    mask[0, : nx // 3] = 0.0  # a stretch of land at the coast

    s_rho = (np.arange(1, N + 1) - N - 0.5) / N
    s_w = (np.arange(N + 1) - N) / N
    v = VERT
    Cs_r = cs(s_rho, v["Vstretching"], v["theta_s"], v["theta_b"])
    Cs_w = cs(s_w, v["Vstretching"], v["theta_s"], v["theta_b"])
    z_rho = depth(s_rho, Cs_r, h, zeta, v["Vtransform"], v["hc"])
    z_w = depth(s_w, Cs_w, h, zeta, v["Vtransform"], v["hc"])
    # a made-up temperature that belongs to the place and the depth
    temp = (2.0 + 16.0 * np.exp(z_rho / 300.0)
            + 2.0 * np.cos(lat * D2R * 6.0) * np.exp(z_rho / 800.0)).astype("float32")
    temp[:, mask == 0] = np.nan

    rho = ("eta_rho", "xi_rho")
    ds = xr.Dataset(
        {
            "temp": (("s_rho",) + rho, temp, {
                "long_name": "potential temperature", "units": "Celsius",
                "coordinates": "lon_rho lat_rho s_rho"}),
            "h": (rho, h, {"long_name": "bathymetry at RHO-points", "units": "meter",
                           "coordinates": "lon_rho lat_rho"}),
            "zeta": (rho, zeta, {"long_name": "free-surface", "units": "meter",
                                 "coordinates": "lon_rho lat_rho"}),
            "mask_rho": (rho, mask, {"long_name": "mask on RHO-points",
                                     "flag_values": [0.0, 1.0], "flag_meanings": "land water",
                                     "coordinates": "lon_rho lat_rho"}),
            "Cs_r": (("s_rho",), Cs_r, {
                "long_name": "S-coordinate stretching curves at RHO-points",
                "valid_min": -1.0, "valid_max": 0.0}),
            "hc": ((), v["hc"], {"long_name": "S-coordinate parameter, critical depth", "units": "meter"}),
            "Vtransform": ((), np.int32(v["Vtransform"]), {"long_name": "vertical terrain-following transformation equation"}),
            "Vstretching": ((), np.int32(v["Vstretching"]), {"long_name": "vertical terrain-following stretching function"}),
            "theta_s": ((), v["theta_s"], {"long_name": "S-coordinate surface control parameter"}),
            "theta_b": ((), v["theta_b"], {"long_name": "S-coordinate bottom control parameter"}),
        },
        coords={
            "lon_rho": (rho, lon, {"long_name": "longitude of RHO-points",
                                   "units": "degree_east", "standard_name": "longitude"}),
            "lat_rho": (rho, lat, {"long_name": "latitude of RHO-points",
                                   "units": "degree_north", "standard_name": "latitude"}),
            "s_rho": (("s_rho",), s_rho, {
                "long_name": "S-coordinate at RHO-points", "valid_min": -1.0, "valid_max": 0.0,
                "positive": "up", "standard_name": "ocean_s_coordinate_g2",
                "formula_terms": "s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc"}),
        },
        attrs={"type": "ROMS/TOMS history file (synthetic)", "Conventions": "CF-1.8",
               "title": "wxylib fixtures I1 roms-fan"},
    )
    enc = {k: {"chunks": ds[k].shape if ds[k].ndim < 2 else ds[k].shape[:-2] + (24, 12)}
           for k in ds.variables if ds[k].ndim >= 1}
    for k in ds.variables:
        # xarray would add its own coordinates attribute; keep ROMS's
        ds[k].encoding.pop("coordinates", None)
    store = os.path.join(OUT, "roms-fan.zarr")
    ds.to_zarr(store, zarr_format=3, mode="w", encoding=enc, consolidated=False)

    name = "roms-fan"
    raw = {}
    for key, a, dt in [("lon_rho", lon, "<f8"), ("lat_rho", lat, "<f8"),
                       ("z_rho", z_rho, "<f8"), ("z_w", z_w, "<f8")]:
        raw[key] = "%s-%s.bin" % (name, key.replace("_", "-"))
        a.astype(dt).tofile(os.path.join(OUT, raw[key]))

    # GDAL: the top layer (k = N - 1) through its geolocation arrays
    top = temp[-1]
    src_tif = os.path.join(OUT, name + "-top.tif")
    lonlat_tif = os.path.join(OUT, name + "-lonlat.tif")
    write_tif(src_tif, [top], gdal.GDT_Float32)
    write_tif(lonlat_tif, [lon, lat], gdal.GDT_Float64)
    vrt = os.path.join(OUT, name + "-top.vrt")
    write_geoloc_vrt(vrt, src_tif, lonlat_tif, ny, nx)
    bounds, size = (150.0, -48.0, 172.0, -23.0), (220, 250)
    gt = [bounds[0], (bounds[2] - bounds[0]) / size[0], 0.0,
          bounds[3], 0.0, -(bounds[3] - bounds[1]) / size[1]]
    results = []
    for m in METHODS:
        out = "/vsimem/%s-%s.tif" % (name, m)
        gdal.Warp(out, vrt, format="GTiff", dstSRS="EPSG:4326", outputBounds=bounds,
                  width=size[0], height=size[1], resampleAlg=m, geoloc=True,
                  errorThreshold=0.125, srcNodata="nan", dstNodata="nan",
                  outputType=gdal.GDT_Float32)
        a = gdal.Open(out).ReadAsArray()
        gdal.Unlink(out)
        f = "%s-top-to-4326-%s.bin" % (name, m)
        a.astype("<f4").tofile(os.path.join(OUT, f))
        results.append({"resampling": m, "file": f,
                        "valid": int(np.isfinite(a).sum())})
    return {
        "name": name, "store": "roms-fan.zarr", "case": "I1",
        "shape": {"s_rho": N, "eta_rho": ny, "xi_rho": nx},
        "construction": "playproj fan(ni=24, lon0=150.4, lat0=-33.2, bearing=105, spread=90, r0=0.9, r1=18) at cell centres, spherical",
        "vertical": dict(VERT, N=N, standard_name="ocean_s_coordinate_g2"),
        "truth": {
            "lon_rho": raw["lon_rho"], "lat_rho": raw["lat_rho"],
            "z_rho": raw["z_rho"], "z_w": raw["z_w"],
            "layout": "little-endian float64, C order; z_rho (s_rho, eta_rho, xi_rho), z_w (N + 1, eta_rho, xi_rho), bottom first as ROMS orders k",
            "z_from": "NumPy port of ROMS set_scoord.F and set_depth.F (playproj test/make-roms-fixture.py), with zeta",
        },
        "warp": {
            "layer": "temp[s_rho = N - 1] (surface)",
            "geolocation": "PIXEL_CENTER", "error_threshold": 0.125,
            "crs": "EPSG:4326", "geotransform": gt, "width": size[0], "height": size[1],
            "dtype": "float32", "nodata": "nan", "results": results,
        },
    }


# ---------------------------------------------------------------- I3

POLE = {"grid_north_pole_latitude": 39.25, "grid_north_pole_longitude": -162.0}


def rotated_lonlat(rlon, rlat):
    """CF rotated pole (PROJ ob_tran) on the sphere, degrees in and out."""
    RL, RP = rlon * D2R, rlat * D2R
    pp = POLE["grid_north_pole_latitude"] * D2R
    lat = np.arcsin(np.cos(pp) * np.cos(RP) * np.cos(RL) + np.sin(pp) * np.sin(RP))
    lon = np.arctan2(np.cos(RP) * np.sin(RL),
                     np.sin(pp) * np.cos(RP) * np.cos(RL) - np.cos(pp) * np.sin(RP))
    return wrap(lon / D2R + POLE["grid_north_pole_longitude"] + 180.0), lat / D2R


def rotated():
    # EURO-CORDEX EUR-11 window, at 1 degree instead of 0.11
    step, rl0, rp0 = 1.0, -28.375, -23.375
    nx, ny = 47, 45
    rlon = rl0 + step * np.arange(nx)
    rlat = rp0 + step * np.arange(ny)
    RL, RP = np.meshgrid(rlon, rlat)
    lon, lat = rotated_lonlat(RL, RP)
    tas = (285.0 - 30.0 * np.sin(lat * D2R) ** 2 + 3.0 * np.cos(3 * lon * D2R)).astype("float32")

    gm = dict(POLE, grid_mapping_name="rotated_latitude_longitude")
    ds = xr.Dataset(
        {
            "tas": (("rlat", "rlon"), tas, {"standard_name": "air_temperature", "units": "K",
                                            "grid_mapping": "rotated_pole",
                                            "coordinates": "lon lat"}),
            "rotated_pole": ((), np.int32(0), gm),
        },
        coords={
            "rlon": ("rlon", rlon, {"standard_name": "grid_longitude", "units": "degrees",
                                    "long_name": "longitude in rotated pole grid", "axis": "X"}),
            "rlat": ("rlat", rlat, {"standard_name": "grid_latitude", "units": "degrees",
                                    "long_name": "latitude in rotated pole grid", "axis": "Y"}),
            "lon": (("rlat", "rlon"), lon.astype("float32"),
                    {"standard_name": "longitude", "units": "degrees_east"}),
            "lat": (("rlat", "rlon"), lat.astype("float32"),
                    {"standard_name": "latitude", "units": "degrees_north"}),
        },
        attrs={"Conventions": "CF-1.8", "title": "wxylib fixtures I3 rotated (EURO-CORDEX window, 1 degree)"},
    )
    for k in ds.variables:
        ds[k].encoding.pop("coordinates", None)
    store = os.path.join(OUT, "rotated.zarr")
    ds.to_zarr(store, zarr_format=3, mode="w", consolidated=False,
               encoding={"tas": {"chunks": (ny, nx)}})

    # the oracle: PROJ's ob_tran (osgeo.osr) on the generating affine
    proj = ("+proj=ob_tran +o_proj=longlat +o_lat_p=%r +o_lon_p=0 +lon_0=%r "
            "+R=6371229 +no_defs" % (POLE["grid_north_pole_latitude"],
                                     POLE["grid_north_pole_longitude"] + 180.0))
    src = osr.SpatialReference()
    src.ImportFromProj4(proj)
    dst = osr.SpatialReference()
    dst.ImportFromProj4("+proj=longlat +R=6371229 +no_defs")
    for s in (src, dst):
        s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    ct = osr.CoordinateTransformation(src, dst)
    pts = np.array(ct.TransformPoints(np.stack([RL.ravel(), RP.ravel()], 1).tolist()))
    plon = pts[:, 0].reshape(RL.shape)
    plat = pts[:, 1].reshape(RL.shape)
    dlon = np.abs(wrap(plon - lon.astype("float32").astype("float64"))).max()
    dlat = np.abs(plat - lat.astype("float32").astype("float64")).max()
    dlon64 = np.abs(wrap(plon - lon)).max()
    dlat64 = np.abs(plat - lat).max()
    assert max(dlon64, dlat64) < 1e-9, (dlon64, dlat64)
    assert max(dlon, dlat) < 1e-5, (dlon, dlat)

    return {
        "name": "rotated", "store": "rotated.zarr", "case": "I3",
        "shape": {"rlat": ny, "rlon": nx},
        "grid_mapping": gm,
        "proj": proj,
        "affine": {
            "frame": "rotated degrees (rlon, rlat), cell centres",
            "rlon": [rl0, step], "rlat": [rp0, step],
            "geotransform": [rl0 - step / 2, step, 0.0, rp0 - step / 2, 0.0, step],
            "note": "GDAL order over (column, row) corners; rlat ascends with the row index",
        },
        "check": {
            "by": "osgeo.osr ob_tran of the affine's centres against the stored 2D lon, lat",
            "max_abs_deg_float32": {"lon": float(dlon), "lat": float(dlat)},
            "max_abs_deg_float64": {"lon": float(dlon64), "lat": float(dlat64)},
        },
    }


def main():
    shutil.rmtree(OUT, ignore_errors=True)
    os.makedirs(OUT)
    cases = [roms_fan(), rotated()]
    man = {"gdal": gdal.VersionInfo("RELEASE_NAME"), "osgeo": osgeo.__version__,
           "zarr_format": 3, "cases": cases}
    with open(os.path.join(OUT, "manifest.json"), "w") as f:
        json.dump(man, f, indent=2)
        f.write("\n")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
