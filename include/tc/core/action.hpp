#pragma once

struct Action {
    int from_idx = -1;
    int gen = -1;

    Action() = default;

    Action(const Action &) = default;

    Action(int from_idx, int gen)
        : from_idx(from_idx), gen(gen) {
    }
};
