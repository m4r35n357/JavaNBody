/**
 * 
 */
package uk.me.doitto.orbits;

/**
 * @author ian
 *
 */
public enum Integrators implements IIntegrators {
	EULER {
		@Override
		public void solve (Symplectic s, Coordinates first, Coordinates second) {
			assert first != second;
			first.update(s, 1.0);
			second.update(s, 1.0);
		}
	},
	STORMER_VERLET_2 {
		@Override
		public void solve (Symplectic s, Coordinates first, Coordinates second) {
			assert first != second;
			first.update(s, 0.5);
			second.update(s, 1.0);
			first.update(s, 0.5);
		}
	},
	STORMER_VERLET_4 {
		@Override
		public void solve (Symplectic s, Coordinates first, Coordinates second) {
			assert first != second;
			double c1, c2, c3, c4, d1, d2, d3;
			double cubeRoot2 = Math.pow(2.0, 1.0 / 3.0);
			double denom = 2.0 - cubeRoot2;
			c1 = c4 = 1.0 / (2.0 * denom);
			c2 = c3 = (1.0 - cubeRoot2) / (2.0 * denom);
			d1 = d3 = 1.0 / denom;
			d2 = - cubeRoot2 / denom;
			first.update(s, c1);
			second.update(s, d1);
			first.update(s, c2);
			second.update(s, d2);
			first.update(s, c3);
			second.update(s, d3);
			first.update(s, c4);
		}
	},
	STORMER_VERLET_8 {
		@Override
		public void solve (Symplectic s, Coordinates first, Coordinates second) {
			assert first != second;
			first.update(s, -0.74167036435061295344822780);
			second.update(s, 0.40910082580003159399730010);
			first.update(s, -0.19075471029623837995387626);
			second.update(s, 0.57386247111608226665638773);
			first.update(s, -0.29906418130365592384446354);
			second.update(s, -0.33462491824529818378495798);
			first.update(s, -0.31529309239676659663205666);
			second.update(s, 0.79688793935291635401978884);
			first.update(s, -0.31529309239676659663205666);
			second.update(s, -0.33462491824529818378495798);
			first.update(s, -0.29906418130365592384446354);
			second.update(s, 0.57386247111608226665638773);
			first.update(s, -0.19075471029623837995387626);
			second.update(s, 0.40910082580003159399730010);
			first.update(s, -0.74167036435061295344822780);
		}
	};

	@Override
	public abstract void solve (Symplectic s, Coordinates first, Coordinates second);
}