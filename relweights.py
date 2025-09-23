"""Utilities for constructing RelWeights matrices and exporting GAL files.

The implementation mirrors the README description: we compute an incidence matrix
between a relational layer (A) and an inherited layer (B), derive the relational
weights R = W W' with a zero diagonal, and export the structure for PySAL workflows.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence

import geopandas as gpd
import numpy as np
import pandas as pd
from libpysal import io as psio
from libpysal import weights


@dataclass
class RelWeightsResult:
    """Container for the matrices and PySAL weight object."""

    rel_matrix: np.ndarray
    weights_matrix: np.ndarray
    w: weights.W


def _prepare_ids(ids: Optional[Sequence], count: int) -> List[str]:
    if ids is None:
        return [str(i) for i in range(count)]
    if len(ids) != count:
        raise ValueError("Length of ids must match the number of rows in rel_layer.")
    return [str(i) for i in ids]


def tab2relweights(
    rel_layer: gpd.GeoDataFrame,
    inh_layer: gpd.GeoDataFrame,
    *,
    binary: bool = True,
    transform: Optional[str] = "R",
    ids: Optional[Iterable] = None,
    output_gal: Optional[str] = None,
) -> RelWeightsResult:
    """Compute RelWeights across two overlay layers and optionally export a GAL file.

    Parameters
    ----------
    rel_layer : GeoDataFrame
        Layer whose features will inherit relational neighbours.
    inh_layer : GeoDataFrame
        Layer that induces inherited contiguities.
    binary : bool, default True
        Convert positive weights to 1 prior to building the PySAL weight object.
    transform : str or None, default "R"
        Row standardisation applied to the PySAL weight via `W.transform`.
    ids : iterable or None
        Optional iterable of feature identifiers aligned with `rel_layer` rows.
    output_gal : str or None
        When provided, the result is exported to a .gal file using libpysal.

    Returns
    -------
    RelWeightsResult
        Dataclass containing the dense relational matrix, the weights matrix actually used
        to create the PySAL `W` object, and the `W` object itself.
    """

    if not isinstance(rel_layer, gpd.GeoDataFrame):
        raise TypeError("rel_layer must be a GeoDataFrame.")
    if not isinstance(inh_layer, gpd.GeoDataFrame):
        raise TypeError("inh_layer must be a GeoDataFrame.")

    rel = rel_layer.reset_index(drop=True).copy()
    inh = inh_layer.reset_index(drop=True).copy()

    rel["id_x"] = np.arange(len(rel))
    inh["id_y"] = np.arange(len(inh))

    intersections = gpd.sjoin(
        rel[["id_x", "geometry"]],
        inh[["id_y", "geometry"]],
        how="inner",
        predicate="intersects",
    )

    rel_ids = _prepare_ids(ids, len(rel))

    if intersections.empty:
        rel_matrix = np.zeros((len(rel), len(rel)), dtype=float)
        weights_matrix = rel_matrix.copy()
    else:
        attrs = intersections[["id_x", "id_y"]]
        incidence = pd.crosstab(attrs["id_x"], attrs["id_y"])
        incidence = incidence.reindex(
            index=np.arange(len(rel)),
            columns=np.arange(len(inh)),
            fill_value=0,
        )
        incidence_values = incidence.to_numpy(dtype=float)
        rel_matrix = incidence_values @ incidence_values.T
        np.fill_diagonal(rel_matrix, 0.0)
        if binary:
            weights_matrix = (rel_matrix > 0).astype(float)
        else:
            weights_matrix = rel_matrix.astype(float)

    w = weights.full2W(weights_matrix, ids=rel_ids)
    if transform:
        w.transform = transform

    if output_gal:
        with psio.open(output_gal, "w") as gal:
            gal.write(w)

    return RelWeightsResult(rel_matrix=rel_matrix, weights_matrix=weights_matrix, w=w)
