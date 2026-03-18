# Boundary Logic — Function Flow Diagram

Paste the code block below into [https://mermaid.live](https://mermaid.live) to render the diagram:

flowchart TD

    %% ─────────────────────────────────────────
    %% INPUT
    %% ─────────────────────────────────────────
    RAW[(Raw data frame\nwith class column)]
    EXT[(External / new\nobservation)]

    %% ─────────────────────────────────────────
    %% PHASE 1 — Setup
    %% ─────────────────────────────────────────
    subgraph P1["PHASE 1 — Model & Projection Setup"]
        direction TB

        A1["bl_prepare_data()\nbl_wrap_data()\n→ bl_data"]
        A2["bl_filter_outliers()\n→ bl_filter_result\n(polygon + retained rows)"]
        A3A["bl_fit_model()\n→ bl_model"]
        A3B["bl_wrap_model()\n→ bl_model\n(custom / GAM / XGB)"]
        A4["bl_build_projection()\n→ bl_projection\n(V, tV, X_center, X_sd,\nbiplot_obj)"]
        A5["bl_build_grid()\n→ bl_grid\n(m×m grid, ct, ct_surrogate,\npolygon, col_value)"]
        A6["bl_assemble()\n→ bl_result\n(single anchor object)"]

        A1 --> A2 --> A3A & A3B
        A3A & A3B --> A4 --> A5 --> A6
        A2 --> A6
    end

    RAW --> A1

    %% ─────────────────────────────────────────
    %% PHASE 1 — Plotting
    %% ─────────────────────────────────────────
    subgraph P1P["PHASE 1 — Visualisation"]
        direction TB

        B1["bl_project_points()\n→ bl_points\n(Z, pred_prob, pred_col)"]
        B2["bl_predict()\n→ data.frame\n(row, pred_prob, pred_class,\ntrue_class, confusion,\nfeatures)"]
        B3["plot_biplotEZ()\n[biplotEZ pipeline:\naxes → grid → points\n→ contours → target\n→ boundary overlay]"]

        B1 --> B3
        B2
    end

    A6 --> B1 & B2 & B3

    %% ─────────────────────────────────────────
    %% PHASE 2 — Global Interpretations
    %% ─────────────────────────────────────────
    subgraph P2["PHASE 2 — Global Interpretations"]
        direction TB

        C1["bl_find_boundary()\n→ bl_boundary\n(Z_obs, B_z, B_x,\ndist_z per observation)\nConstraints:\n• hull polygon\n• train_ranges"]
        C2["plot.bl_boundary()\njitter / boxplot\nTP/TN/FP/FN colours"]
        C3["bl_robustness()\n→ per-variable distance\n+ scalar robustness score"]
        C4["bl_surrogate()\n→ bl_surrogate\n(region assignment,\naccuracy vs model\n+ vs true labels)"]
        C5["plot.bl_surrogate()\nhull-clipped grid\n+ surrogate-coloured points"]

        C1 --> C2 & C3
        C4 --> C5
    end

    A6 --> C1 & C4
    B1 -.->|"points = bl_project_points(test_data)"| B3

    %% ─────────────────────────────────────────
    %% PHASE 3 — Local Interpretations
    %% ─────────────────────────────────────────
    subgraph P3["PHASE 3 — Local Interpretations"]
        direction TB

        D1["bl_select_target()\n→ bl_target\n(x_obs, z_obs,\npred_prob, pred_class)"]
        D2["set_filters()\n→ bl_filters\ndecrease / increase /\nfixed ±0.5 / c(min,max)"]
        D3["bl_find_boundary_local()\n→ bl_local_result\nFor each eigenvector pair:\n• SVD rotation (.bl_rotate)\n• local m×m grid\n• contour extraction\n• train_ranges filter\n• actionability filter\nPicks pair with min dist_z\nConstraints:\n• train_ranges\n• set_filters\n• NO hull polygon"]
        D4["plot.bl_local_result()\n[biplotEZ pipeline\npatched with Vr_rot /\nZ_train_rot / ax.one.unit:\naxes → grid → points\n→ contours → target\n→ CF cross + arrow]"]
        D5["bl_shapley()\n→ bl_shapley\nExact (p≤14) or\npermutation Shapley\nSupports / Contradicts\nper variable"]
        D6["plot.bl_shapley()\nggplot2 bar chart\ndarkblue=Supports\ngrey=Contradicts"]
        D7["bl_sparse_cf()\n→ bl_sparse_result\n• keep Supports at CF value\n• revert Contradicts to obs\n• revert fixed to obs\n• round to nearest round_to\n• re-score → solution_valid"]
        D8["plot.bl_sparse_result()\nbiplotEZ local biplot +\ngrey cross = full CF\ngreen = sparse CF valid\nyellow = sparse CF invalid"]
        D9["print.bl_sparse_result()\nObserved pred\nFull CF pred\nSparse CF pred"]

        D1 --> D2 --> D3
        D3 --> D4
        D3 --> D5 --> D6
        D5 --> D7 --> D8 & D9
    end

    A6 --> D1
    EXT -->|"single-row\ndata frame"| D1

    %% ─────────────────────────────────────────
    %% Cross-phase links
    %% ─────────────────────────────────────────
    C1 -.->|"boundary = bl_bnd"| B3
    B2 -.->|"inspect → choose tdp"| D1

    %% ─────────────────────────────────────────
    %% Styling
    %% ─────────────────────────────────────────
    classDef phase1    fill:#dbeafe,stroke:#3b82f6,color:#1e3a5f
    classDef phase1p   fill:#eff6ff,stroke:#93c5fd,color:#1e3a5f
    classDef phase2    fill:#dcfce7,stroke:#22c55e,color:#14532d
    classDef phase3    fill:#fef9c3,stroke:#eab308,color:#713f12
    classDef data      fill:#f3e8ff,stroke:#a855f7,color:#3b0764
    classDef external  fill:#ffe4e6,stroke:#f43f5e,color:#881337

    class A1,A2,A3A,A3B,A4,A5,A6 phase1
    class B1,B2,B3 phase1p
    class C1,C2,C3,C4,C5 phase2
    class D1,D2,D3,D4,D5,D6,D7,D8,D9 phase3
    class RAW data
    class EXT external

