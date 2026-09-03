#include "dmft_loop_bridge.hpp"

#include "../legacy/headers.new.h"
#include "../legacy/parameters.cpp"
#include "../legacy/solver.cpp"

namespace legacy_bridge {

DmftLoopSnapshot run(double lambda0, double hz, double factor, unsigned int seed) {
    INT_RANDOM = static_cast<int>(seed);

    DMFT_AF_solver solver(lambda0);
    solver.init(hz);

    DmftLoopSnapshot snap;
    snap.raw_loops_result = solver.raw_loops(factor);
    snap.fine_loops_result = solver.fine_loops(factor);

    snap.gu.assign(solver.Sf.gu, solver.Sf.gu + Nw);
    snap.gd.assign(solver.Sf.gd, solver.Sf.gd + Nw);
    snap.delta_u.assign(solver.Sf.Delta_u, solver.Sf.Delta_u + Nw);
    snap.delta_d.assign(solver.Sf.Delta_d, solver.Sf.Delta_d + Nw);
    snap.hz_fine = solver.Sf.hz_fine;
    return snap;
}

}  // namespace legacy_bridge
