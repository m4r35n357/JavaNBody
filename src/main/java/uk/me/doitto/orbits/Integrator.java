package uk.me.doitto.orbits;

/**
 * @author ian
 * <p>
 * Symplectic integrators of various orders and types
 */
public enum Integrator {
	/**
	 * Euler 1st-order
	 */
	EULER {
		@Override
		void solve (Symplectic s, PhaseSpace first, PhaseSpace second) {
			first.update(s, 1.0);
			second.update(s, 1.0);
		}
	},
	/**
	 * Stormer-Verlet 2nd-order
	 */
	STORMER_VERLET_2 {
		@Override
		void solve (Symplectic s, PhaseSpace first, PhaseSpace second) {
			symplectic2 (s, first, second, 1.0);
		}
	},
	/**
	 * Yoshida 4th-order
	 */
	STORMER_VERLET_4 {
		@Override
		void solve (Symplectic s, PhaseSpace first, PhaseSpace second) {
			symplectic2(s, first, second, 1.351207191959657);
			symplectic2(s, first, second, -1.702414383919315);
			symplectic2(s, first, second, 1.351207191959657);
		}
	},
	/**
	 * Yoshida 6th-order
	 */
	STORMER_VERLET_6 {
		@Override
		void solve (Symplectic s, PhaseSpace first, PhaseSpace second) {
			symplectic2(s, first, second, 0.784513610477560e0);
			symplectic2(s, first, second, 0.235573213359357e0);
			symplectic2(s, first, second, -1.17767998417887e0);
			symplectic2(s, first, second, 1.31518632068391e0);
			symplectic2(s, first, second, -1.17767998417887e0);
			symplectic2(s, first, second, 0.235573213359357e0);
			symplectic2(s, first, second, 0.784513610477560e0);
		}
	},
	/**
	 * Yoshida 8th-order
	 */
	STORMER_VERLET_8 {
		@Override
		void solve (Symplectic s, PhaseSpace first, PhaseSpace second) {
			symplectic2(s, first, second, 0.104242620869991e1);
			symplectic2(s, first, second, 0.182020630970714e1);
			symplectic2(s, first, second, 0.157739928123617e0);
			symplectic2(s, first, second, 0.244002732616735e1);
			symplectic2(s, first, second, -0.716989419708120e-2);
			symplectic2(s, first, second, -0.244699182370524e1);
			symplectic2(s, first, second, -0.161582374150097e1);
			symplectic2(s, first, second, -0.17808286265894516e1);
			symplectic2(s, first, second, -0.161582374150097e1);
			symplectic2(s, first, second, -0.244699182370524e1);
			symplectic2(s, first, second, -0.716989419708120e-2);
			symplectic2(s, first, second, 0.244002732616735e1);
			symplectic2(s, first, second, 0.157739928123617e0);
			symplectic2(s, first, second, 0.182020630970714e1);
			symplectic2(s, first, second, 0.104242620869991e1);
		}
	};

	/**
	 * Basic 2nd-order Stormer-Verlet step which is composed into higher order methods
	 */
	protected final void symplectic2 (Symplectic s, PhaseSpace first, PhaseSpace second, double step) {
		first.update(s, 0.5 * step);
		second.update(s, step);
		first.update(s, 0.5 * step);
	}
	
	/**
	 * Perform one iteration step for the configured integrator
	 * @param s the Symplectic instance
	 * @param first the "outer" update
	 * @param second the "middle" update
	 */
	abstract void solve (Symplectic s, PhaseSpace first, PhaseSpace second);
}
