// Copyright (C) 2026 Shingo Hirano and Sho Higashi
// Licensed under the MIT found in the
// https://github.com/astro-sim-lab/arche/blob/main/LICENSE
#pragma once
// ---------------------------------------------------------------------------
// model_traits.h — aggregator of the per-model compile-time trait structs
//
// Each model's trait lives next to its definition in models/<model>/traits.h.
// This header bundles them so kernel templates that are generic over the model
// can include a single trait entry point.
//
//   Nakauchi2019, Nakauchi2019_Minimal  → models/primordial/traits.h
//   Nakauchi2021, Nakauchi2021_Minimal  → models/metal_grain/traits.h
// ---------------------------------------------------------------------------
#include "models/metal_grain/traits.h"
#include "models/primordial/traits.h"
