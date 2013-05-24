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
			sympBase (s, first, second, 1.0);
		}
	},
	/**
	 * Yoshida 4th-order
	 */
	STORMER_VERLET_4 {
		@Override
		void solve (Symplectic s, PhaseSpace first, PhaseSpace second) {
			double y = 1.0 / (2.0 - CUBE_ROOT_2);
			sympBase(s, first, second, y);
			sympBase(s, first, second, - CUBE_ROOT_2 * y);
			sympBase(s, first, second, y);
		}
	},
	/**
	 * Yoshida 6th-order
	 */
	STORMER_VERLET_6 {
		@Override
		void solve (Symplectic s, PhaseSpace first, PhaseSpace second) {
			sympBase(s, first, second, 0.78451361047755726381949763);
			sympBase(s, first, second, 0.23557321335935813368479318);
			sympBase(s, first, second, -1.17767998417887100694641568);
			sympBase(s, first, second, 1.31518632068391121888424973);
			sympBase(s, first, second, -1.17767998417887100694641568);
			sympBase(s, first, second, 0.23557321335935813368479318);
			sympBase(s, first, second, 0.78451361047755726381949763);
		}
	},
	/**
	 * Yoshida 8th-order
	 */
	STORMER_VERLET_8 {
		@Override
		void solve (Symplectic s, PhaseSpace first, PhaseSpace second) {
			sympBase(s, first, second, 0.74167036435061295344822780);
			sympBase(s, first, second, -0.40910082580003159399730010);
			sympBase(s, first, second, 0.19075471029623837995387626);
			sympBase(s, first, second, -0.57386247111608226665638773);
			sympBase(s, first, second, 0.29906418130365592384446354);
			sympBase(s, first, second, 0.33462491824529818378495798);
			sympBase(s, first, second, 0.31529309239676659663205666);
			sympBase(s, first, second, -0.79688793935291635401978884);
			sympBase(s, first, second, 0.31529309239676659663205666);
			sympBase(s, first, second, 0.33462491824529818378495798);
			sympBase(s, first, second, 0.29906418130365592384446354);
			sympBase(s, first, second, -0.57386247111608226665638773);
			sympBase(s, first, second, 0.19075471029623837995387626);
			sympBase(s, first, second, -0.40910082580003159399730010);
			sympBase(s, first, second, 0.74167036435061295344822780);
		}
	},
	/**
	 * Yoshida 10th-order
	 */
	STORMER_VERLET_10 {
		@Override
		void solve (Symplectic s, PhaseSpace first, PhaseSpace second) {
			sympBase(s, first, second, 0.09040619368607278492161150);
			sympBase(s, first, second, 0.53591815953030120213784983);
			sympBase(s, first, second, 0.35123257547493978187517736);
			sympBase(s, first, second, -0.31116802097815835426086544);
			sympBase(s, first, second, -0.52556314194263510431065549);
			sympBase(s, first, second, 0.14447909410225247647345695);
			sympBase(s, first, second, 0.02983588609748235818064083);
			sympBase(s, first, second, 0.17786179923739805133592238);
			sympBase(s, first, second, 0.09826906939341637652532377);
			sympBase(s, first, second, 0.46179986210411860873242126);
			sympBase(s, first, second, -0.33377845599881851314531820);
			sympBase(s, first, second, 0.07095684836524793621031152);
			sympBase(s, first, second, 0.23666960070126868771909819);
			sympBase(s, first, second, -0.49725977950660985445028388);
			sympBase(s, first, second, -0.30399616617237257346546356);
			sympBase(s, first, second, 0.05246957188100069574521612);
			sympBase(s, first, second, 0.44373380805019087955111365);
			sympBase(s, first, second, 0.05246957188100069574521612);
			sympBase(s, first, second, -0.30399616617237257346546356);
			sympBase(s, first, second, -0.49725977950660985445028388);
			sympBase(s, first, second, 0.23666960070126868771909819);
			sympBase(s, first, second, 0.07095684836524793621031152);
			sympBase(s, first, second, -0.33377845599881851314531820);
			sympBase(s, first, second, 0.46179986210411860873242126);
			sympBase(s, first, second, 0.09826906939341637652532377);
			sympBase(s, first, second, 0.17786179923739805133592238);
			sympBase(s, first, second, 0.02983588609748235818064083);
			sympBase(s, first, second, 0.14447909410225247647345695);
			sympBase(s, first, second, -0.52556314194263510431065549);
			sympBase(s, first, second, -0.31116802097815835426086544);
			sympBase(s, first, second, 0.35123257547493978187517736);
			sympBase(s, first, second, 0.53591815953030120213784983);
			sympBase(s, first, second, 0.09040619368607278492161150);
		}
	};

	private static final double CUBE_ROOT_2 = Math.pow(2.0, 1.0 / 3.0);
	
	/**
	 * Basic 2nd-order Stormer-Verlet step which is composed into higher order methods
	 */
	protected final void sympBase (Symplectic s, PhaseSpace first, PhaseSpace second, double c) {
		first.update(s, 0.5 * c);
		second.update(s, c);
		first.update(s, 0.5 * c);
	}
	
	/**
	 * Perform one iteration step for the configured integrator
	 * @param s the Symplectic instance
	 * @param first the "outer" update
	 * @param second the "middle" update
	 */
	abstract void solve (Symplectic s, PhaseSpace first, PhaseSpace second);
}
