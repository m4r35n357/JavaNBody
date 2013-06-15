package uk.me.doitto.orbits;

import java.math.BigDecimal;

/**
 * @author ian
 * <p>
 * Symplectic integrators of various orders and types
 */
public enum Integrator {
	/**
	 * Stormer-Verlet 2nd-order
	 */
	STORMER_VERLET_2 {
		@Override
		void init() {
			BigDecimal[] coefficients = new BigDecimal[1];
			coefficients[0] = new BigDecimal("1.0", Symplectic.P);
			this.coefficients = coefficients;
		}
	},
	/**
	 * Yoshida 4th-order
	 */
	STORMER_VERLET_4 {
		@Override
		void init() {
			BigDecimal CUBE_ROOT_2 = BigDecimal.valueOf(Math.pow(2.0, 1.0 / 3.0));
			BigDecimal y = new BigDecimal("1.0", Symplectic.P).divide(new BigDecimal("2.0", Symplectic.P).subtract(CUBE_ROOT_2), Symplectic.P);
			BigDecimal[] coefficients = new BigDecimal[2];
			coefficients[0] = y;
			coefficients[1] = y.multiply(CUBE_ROOT_2, Symplectic.P).negate(Symplectic.P);
			this.coefficients = coefficients;
		}
	},
	/**
	 * Yoshida 6th-order
	 */
	STORMER_VERLET_6 {
		@Override
		void init() {
			BigDecimal[] coefficients = new BigDecimal[4];
			coefficients[0] = new BigDecimal("0.78451361047755726381949763", Symplectic.P);
			coefficients[1] = new BigDecimal("0.23557321335935813368479318", Symplectic.P);
			coefficients[2] = new BigDecimal("-1.17767998417887100694641568", Symplectic.P);
			coefficients[3] = new BigDecimal("1.31518632068391121888424973", Symplectic.P);
			this.coefficients = coefficients;
		}
	},
	/**
	 * Yoshida 8th-order
	 */
	STORMER_VERLET_8 {
		@Override
		void init() {
			BigDecimal[] coefficients = new BigDecimal[8];
			coefficients[0] = new BigDecimal("0.74167036435061295344822780", Symplectic.P);
			coefficients[1] = new BigDecimal("-0.40910082580003159399730010", Symplectic.P);
			coefficients[2] = new BigDecimal("0.19075471029623837995387626", Symplectic.P);
			coefficients[3] = new BigDecimal("-0.57386247111608226665638773", Symplectic.P);
			coefficients[4] = new BigDecimal("0.29906418130365592384446354", Symplectic.P);
			coefficients[5] = new BigDecimal("0.33462491824529818378495798", Symplectic.P);
			coefficients[6] = new BigDecimal("0.31529309239676659663205666", Symplectic.P);
			coefficients[7] = new BigDecimal("-0.79688793935291635401978884", Symplectic.P);
			this.coefficients = coefficients;
		}
	},
	/**
	 * Yoshida 10th-order
	 */
	STORMER_VERLET_10 {
		@Override
		void init() {
			BigDecimal[] coefficients = new BigDecimal[17];
			coefficients[0] = new BigDecimal("0.09040619368607278492161150", Symplectic.P);
			coefficients[1] = new BigDecimal("0.53591815953030120213784983", Symplectic.P);
			coefficients[2] = new BigDecimal("0.35123257547493978187517736", Symplectic.P);
			coefficients[3] = new BigDecimal("-0.31116802097815835426086544", Symplectic.P);
			coefficients[4] = new BigDecimal("-0.52556314194263510431065549", Symplectic.P);
			coefficients[5] = new BigDecimal("0.14447909410225247647345695", Symplectic.P);
			coefficients[6] = new BigDecimal("0.02983588609748235818064083", Symplectic.P);
			coefficients[7] = new BigDecimal("0.17786179923739805133592238", Symplectic.P);
			coefficients[8] = new BigDecimal("0.09826906939341637652532377", Symplectic.P);
			coefficients[9] = new BigDecimal("0.46179986210411860873242126", Symplectic.P);
			coefficients[10] = new BigDecimal("-0.33377845599881851314531820", Symplectic.P);
			coefficients[11] = new BigDecimal("0.07095684836524793621031152", Symplectic.P);
			coefficients[12] = new BigDecimal("0.23666960070126868771909819", Symplectic.P);
			coefficients[13] = new BigDecimal("-0.49725977950660985445028388", Symplectic.P);
			coefficients[14] = new BigDecimal("-0.30399616617237257346546356", Symplectic.P);
			coefficients[15] = new BigDecimal("0.05246957188100069574521612", Symplectic.P);
			coefficients[16] = new BigDecimal("0.44373380805019087955111365", Symplectic.P);
			this.coefficients = coefficients;
		}
	};

	protected BigDecimal[] coefficients;
	
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
	 * Set up the integrator composition coefficients
	 */
	abstract void init ();
	
	/**
	 * Perform one iteration step for the configured integrator
	 * @param s the Symplectic object reference, for passing through to the Q & P update methods
	 */
	void solve (Symplectic s) {
		int tmp = coefficients.length - 1;
		for (int i = 0; i < tmp; i++) {
			sympBase(s, coefficients[i]);
		}
		for (int i = tmp; i >= 0; i--) {
			sympBase(s, coefficients[i]);
		}
	}
}
