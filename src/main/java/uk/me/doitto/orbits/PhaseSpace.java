package uk.me.doitto.orbits;

/**
 * @author ian
 * 
 * Definitions of the two Hamiltonian coordinates, position and momentum, specifically their update methods
 */
public enum PhaseSpace {	
	Q {  // Position
		@Override
		void update (Symplectic s, double c) {
			for (int i = 0; i < s.np; i++) {
				Particle a = s.particles.get(i);
				double tmp = c / a.mass * s.timeStep;
				a.qX += a.pX * tmp;
				a.qY += a.pY * tmp;
				a.qZ += a.pZ * tmp;
			}
		}
	},	
	P {  // Momentum
		@Override
		void update (Symplectic s, double c) {
			for (int i = 0; i < s.np; i++) {
				Particle a = s.particles.get(i);
				for (int j = 0; j < s.np; j++) {
					Particle b = s.particles.get(j);
					if (i > j) {
						double tmp = - c * s.g * a.mass * b.mass / Math.pow(s.distance(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ), 3) * s.timeStep;
						double dPx = (b.qX - a.qX) * tmp;
						double dPy = (b.qY - a.qY) * tmp;
						double dPz = (b.qZ - a.qZ) * tmp;
						a.pX -= dPx;
						a.pY -= dPy;
						a.pZ -= dPz;
						b.pX += dPx;
						b.pY += dPy;
						b.pZ += dPz;
					}
				}
			}
		}
	};

	abstract void update (Symplectic symplectic, double coefficient);
}
