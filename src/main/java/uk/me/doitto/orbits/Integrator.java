package uk.me.doitto.orbits;

import java.math.BigDecimal;

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
		void solve (Symplectic s) {
			s.updateQ(BIG_ONE);
			s.updateP(BIG_ONE);
		}
	},
	/**
	 * Stormer-Verlet 2nd-order
	 */
	STORMER_VERLET_2 {
		@Override
		void solve (Symplectic s) {
			sympBase (s, new BigDecimal("1.0", Symplectic.P));
		}
	},
	/**
	 * Yoshida 4th-order
	 */
	STORMER_VERLET_4 {
		@Override
		void solve (Symplectic s) {
			BigDecimal y = BIG_ONE.divide(BIG_TW0.subtract(CUBE_ROOT_2), Symplectic.P);
			sympBase(s, y);
			sympBase(s, y.multiply(CUBE_ROOT_2, Symplectic.P).negate(Symplectic.P));
			sympBase(s, y);
		}
	},
	/**
	 * Yoshida 6th-order
	 */
	STORMER_VERLET_6 {
		@Override
		void solve (Symplectic s) {
			sympBase(s, new BigDecimal("0.78451361047755726381949763", Symplectic.P));
			sympBase(s, new BigDecimal("0.23557321335935813368479318", Symplectic.P));
			sympBase(s, new BigDecimal("-1.17767998417887100694641568", Symplectic.P));
			sympBase(s, new BigDecimal("1.31518632068391121888424973", Symplectic.P));
			sympBase(s, new BigDecimal("-1.17767998417887100694641568", Symplectic.P));
			sympBase(s, new BigDecimal("0.23557321335935813368479318", Symplectic.P));
			sympBase(s, new BigDecimal("0.78451361047755726381949763", Symplectic.P));
		}
	},
	/**
	 * Yoshida 8th-order
	 */
	STORMER_VERLET_8 {
		@Override
		void solve (Symplectic s) {
			sympBase(s, new BigDecimal("0.74167036435061295344822780", Symplectic.P));
			sympBase(s, new BigDecimal("-0.40910082580003159399730010", Symplectic.P));
			sympBase(s, new BigDecimal("0.19075471029623837995387626", Symplectic.P));
			sympBase(s, new BigDecimal("-0.57386247111608226665638773", Symplectic.P));
			sympBase(s, new BigDecimal("0.29906418130365592384446354", Symplectic.P));
			sympBase(s, new BigDecimal("0.33462491824529818378495798", Symplectic.P));
			sympBase(s, new BigDecimal("0.31529309239676659663205666", Symplectic.P));
			sympBase(s, new BigDecimal("-0.79688793935291635401978884", Symplectic.P));
			sympBase(s, new BigDecimal("0.31529309239676659663205666", Symplectic.P));
			sympBase(s, new BigDecimal("0.33462491824529818378495798", Symplectic.P));
			sympBase(s, new BigDecimal("0.29906418130365592384446354", Symplectic.P));
			sympBase(s, new BigDecimal("-0.57386247111608226665638773", Symplectic.P));
			sympBase(s, new BigDecimal("0.19075471029623837995387626", Symplectic.P));
			sympBase(s, new BigDecimal("-0.40910082580003159399730010", Symplectic.P));
			sympBase(s, new BigDecimal("0.74167036435061295344822780", Symplectic.P));
		}
	},
	/**
	 * Yoshida 10th-order
	 */
	STORMER_VERLET_10 {
		@Override
		void solve (Symplectic s) {
			sympBase(s, new BigDecimal("0.09040619368607278492161150", Symplectic.P));
			sympBase(s, new BigDecimal("0.53591815953030120213784983", Symplectic.P));
			sympBase(s, new BigDecimal("0.35123257547493978187517736", Symplectic.P));
			sympBase(s, new BigDecimal("-0.31116802097815835426086544", Symplectic.P));
			sympBase(s, new BigDecimal("-0.52556314194263510431065549", Symplectic.P));
			sympBase(s, new BigDecimal("0.14447909410225247647345695", Symplectic.P));
			sympBase(s, new BigDecimal("0.02983588609748235818064083", Symplectic.P));
			sympBase(s, new BigDecimal("0.17786179923739805133592238", Symplectic.P));
			sympBase(s, new BigDecimal("0.09826906939341637652532377", Symplectic.P));
			sympBase(s, new BigDecimal("0.46179986210411860873242126", Symplectic.P));
			sympBase(s, new BigDecimal("-0.33377845599881851314531820", Symplectic.P));
			sympBase(s, new BigDecimal("0.07095684836524793621031152", Symplectic.P));
			sympBase(s, new BigDecimal("0.23666960070126868771909819", Symplectic.P));
			sympBase(s, new BigDecimal("-0.49725977950660985445028388", Symplectic.P));
			sympBase(s, new BigDecimal("-0.30399616617237257346546356", Symplectic.P));
			sympBase(s, new BigDecimal("0.05246957188100069574521612", Symplectic.P));
			sympBase(s, new BigDecimal("0.44373380805019087955111365", Symplectic.P));
			sympBase(s, new BigDecimal("0.05246957188100069574521612", Symplectic.P));
			sympBase(s, new BigDecimal("-0.30399616617237257346546356", Symplectic.P));
			sympBase(s, new BigDecimal("-0.49725977950660985445028388", Symplectic.P));
			sympBase(s, new BigDecimal("0.23666960070126868771909819", Symplectic.P));
			sympBase(s, new BigDecimal("0.07095684836524793621031152", Symplectic.P));
			sympBase(s, new BigDecimal("-0.33377845599881851314531820", Symplectic.P));
			sympBase(s, new BigDecimal("0.46179986210411860873242126", Symplectic.P));
			sympBase(s, new BigDecimal("0.09826906939341637652532377", Symplectic.P));
			sympBase(s, new BigDecimal("0.17786179923739805133592238", Symplectic.P));
			sympBase(s, new BigDecimal("0.02983588609748235818064083", Symplectic.P));
			sympBase(s, new BigDecimal("0.14447909410225247647345695", Symplectic.P));
			sympBase(s, new BigDecimal("-0.52556314194263510431065549", Symplectic.P));
			sympBase(s, new BigDecimal("-0.31116802097815835426086544", Symplectic.P));
			sympBase(s, new BigDecimal("0.35123257547493978187517736", Symplectic.P));
			sympBase(s, new BigDecimal("0.53591815953030120213784983", Symplectic.P));
			sympBase(s, new BigDecimal("0.09040619368607278492161150", Symplectic.P));
		}
	};

	private static final BigDecimal BIG_ONE = new BigDecimal("1.0", Symplectic.P);

	private static final BigDecimal BIG_TW0 = new BigDecimal("2.0", Symplectic.P);

	private static final BigDecimal CUBE_ROOT_2 = BigDecimal.valueOf(Math.pow(2.0, 1.0 / 3.0));
	
	/**
	 * Basic 2nd-order Stormer-Verlet step which is composed into higher order methods
	 * @param s the Symplectic instance
	 * @param c composition coefficient
	 */
	protected final void sympBase (Symplectic s, BigDecimal y) {
		BigDecimal halfY = new BigDecimal("0.5", Symplectic.P).multiply(y, Symplectic.P);
		s.updateQ(halfY);
		s.updateP(y);
		s.updateQ(halfY);
	}
	
	/**
	 * Perform one iteration step for the configured integrator
	 * @param s the Symplectic object reference, for passing through to the Q & P update methods
	 */
	abstract void solve (Symplectic s);
}
