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
			double qr2 = Math.pow(2.0, 1.0 / 3.0);
			double gamma1 = 1.0 / (2.0 - qr2);
			symplectic2(s, first, second, gamma1);
			symplectic2(s, first, second, -qr2 * gamma1);
			symplectic2(s, first, second, gamma1);
		}
	},
	/**
	 * Yoshida 6th-order
	 */
	STORMER_VERLET_6 {
		@Override
		void solve (Symplectic s, PhaseSpace first, PhaseSpace second) {
			symplectic2(s, first, second, 0.78451361047755726381949763);
			symplectic2(s, first, second, 0.23557321335935813368479318);
			symplectic2(s, first, second, -1.17767998417887100694641568);
			symplectic2(s, first, second, 1.31518632068391121888424973);
			symplectic2(s, first, second, -1.17767998417887100694641568);
			symplectic2(s, first, second, 0.23557321335935813368479318);
			symplectic2(s, first, second, 0.78451361047755726381949763);
		}
	},
	/**
	 * Yoshida 8th-order
	 */
	STORMER_VERLET_8 {
		@Override
		void solve (Symplectic s, PhaseSpace first, PhaseSpace second) {
			symplectic2(s, first, second, 0.74167036435061295344822780);
			symplectic2(s, first, second, -0.40910082580003159399730010);
			symplectic2(s, first, second, 0.19075471029623837995387626);
			symplectic2(s, first, second, -0.57386247111608226665638773);
			symplectic2(s, first, second, 0.29906418130365592384446354);
			symplectic2(s, first, second, 0.33462491824529818378495798);
			symplectic2(s, first, second, 0.31529309239676659663205666);
			symplectic2(s, first, second, -0.79688793935291635401978884);
			symplectic2(s, first, second, 0.31529309239676659663205666);
			symplectic2(s, first, second, 0.33462491824529818378495798);
			symplectic2(s, first, second, 0.29906418130365592384446354);
			symplectic2(s, first, second, -0.57386247111608226665638773);
			symplectic2(s, first, second, 0.19075471029623837995387626);
			symplectic2(s, first, second, -0.40910082580003159399730010);
			symplectic2(s, first, second, 0.74167036435061295344822780);
		}
	},
	/**
	 * Yoshida 10th-order
	 */
	STORMER_VERLET_10 {
		@Override
		void solve (Symplectic s, PhaseSpace first, PhaseSpace second) {
			symplectic2(s, first, second, 0.09040619368607278492161150);
			symplectic2(s, first, second, 0.53591815953030120213784983);
			symplectic2(s, first, second, 0.35123257547493978187517736);
			symplectic2(s, first, second, -0.31116802097815835426086544);
			symplectic2(s, first, second, -0.52556314194263510431065549);
			symplectic2(s, first, second, 0.14447909410225247647345695);
			symplectic2(s, first, second, 0.02983588609748235818064083);
			symplectic2(s, first, second, 0.17786179923739805133592238);
			symplectic2(s, first, second, 0.09826906939341637652532377);
			symplectic2(s, first, second, 0.46179986210411860873242126);
			symplectic2(s, first, second, -0.33377845599881851314531820);
			symplectic2(s, first, second, 0.07095684836524793621031152);
			symplectic2(s, first, second, 0.23666960070126868771909819);
			symplectic2(s, first, second, -0.49725977950660985445028388);
			symplectic2(s, first, second, -0.30399616617237257346546356);
			symplectic2(s, first, second, 0.05246957188100069574521612);
			symplectic2(s, first, second, 0.44373380805019087955111365);
			symplectic2(s, first, second, 0.05246957188100069574521612);
			symplectic2(s, first, second, -0.30399616617237257346546356);
			symplectic2(s, first, second, -0.49725977950660985445028388);
			symplectic2(s, first, second, 0.23666960070126868771909819);
			symplectic2(s, first, second, 0.07095684836524793621031152);
			symplectic2(s, first, second, -0.33377845599881851314531820);
			symplectic2(s, first, second, 0.46179986210411860873242126);
			symplectic2(s, first, second, 0.09826906939341637652532377);
			symplectic2(s, first, second, 0.17786179923739805133592238);
			symplectic2(s, first, second, 0.02983588609748235818064083);
			symplectic2(s, first, second, 0.14447909410225247647345695);
			symplectic2(s, first, second, -0.52556314194263510431065549);
			symplectic2(s, first, second, -0.31116802097815835426086544);
			symplectic2(s, first, second, 0.35123257547493978187517736);
			symplectic2(s, first, second, 0.53591815953030120213784983);
			symplectic2(s, first, second, 0.09040619368607278492161150);
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
