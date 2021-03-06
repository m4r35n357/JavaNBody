package uk.me.doitto.orbits;

import static uk.me.doitto.orbits.Integrator.STORMER_VERLET_10;
import static uk.me.doitto.orbits.Integrator.STORMER_VERLET_2;
import static uk.me.doitto.orbits.Integrator.STORMER_VERLET_4;
import static uk.me.doitto.orbits.Integrator.STORMER_VERLET_6;
import static uk.me.doitto.orbits.Integrator.STORMER_VERLET_8;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.json.simple.JSONObject;
import org.json.simple.JSONValue;

/**
 * @author ian
 * <p>
 * Top level class for symplectic integrator simulations
 */
public class Symplectic implements ISymplectic {
	
	public final double iterations, g, timeStep, errorLimit;
	
	public final int np;
	
	public final List<Particle> particles;
	
	private final Integrator integrator;
	
	/**
	 * For calling from JSON data
	 */
	public Symplectic (double g, double simulationTime, double timeStep, double errorLimit, List<Particle> bodies, int integratorOrder) {
		this.particles = new ArrayList<Particle>(bodies);
		this.np = bodies.size();
		this.g = g;
		this.timeStep = timeStep;
		this.errorLimit = errorLimit;
		this.iterations = simulationTime / timeStep;
		switch (integratorOrder) {
		case 2:
			integrator = STORMER_VERLET_2;
			break;
		case 4:
			integrator = STORMER_VERLET_4;
			break;
		case 6:
			integrator = STORMER_VERLET_6;
			break;
		case 8:
			integrator = STORMER_VERLET_8;
			break;
		case 10:
			integrator = STORMER_VERLET_10;
			break;
		default:
			integrator = STORMER_VERLET_4;
			break;
		}
		integrator.init();
	}

	/**
	 * Separation between points A and B
	 * @param xA x coordinate of point A
	 * @param yA y coordinate of point A
	 * @param zA z coordinate of point A
	 * @param xB x coordinate of point B
	 * @param yB y coordinate of point B
	 * @param zB z coordinate of point B
	 * @return the Euclidean distance between points A and B
	 */
	double distance (double xA, double yA, double zA, double xB, double yB, double zB) {
		return Math.sqrt(Math.pow(xB - xA, 2) + Math.pow(yB - yA, 2) + Math.pow(zB - zA, 2));
	}
	
	/**
	 * Total (kinetic + potential) energy of the system
	 * @return the total energy
	 */
	public double hamiltonian () {
		double energy = 0.0;
		for (int i = 0; i < np; i++) {
			Particle a = particles.get(i);
			energy += (0.5 * (a.pX * a.pX + a.pY * a.pY + a.pZ * a.pZ)) / a.mass;
			for (int j = 0; j < np; j++) {
				if (i > j) {
					Particle b = particles.get(j);
					energy -= (g * a.mass * b.mass) / (distance(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ));
				}
			}
		}
		return energy;
	}
	
	/**
	 * Position update implements dH/dp, which in this case is a function of p only
	 * @param c composition coefficient
	 */
	public void updateQ (double c) {
		for (int i = 0; i < np; i++) {
			Particle a = particles.get(i);
			double tmp = c * timeStep / a.mass;
			a.qX += a.pX * tmp;
			a.qY += a.pY * tmp;
			a.qZ += a.pZ * tmp;
		}
	}
	
	/**
	 * Momentum update implements -dH/dq, which in this case is a function of q only
	 * @param c composition coefficient
	 */
	public void updateP (double c) {
		for (int i = 0; i < np; i++) {
			Particle a = particles.get(i);
			for (int j = 0; j < np; j++) {
				if (i > j) {
					Particle b = particles.get(j);
					double tmp = - c * g * a.mass * b.mass * timeStep / Math.pow(distance(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ), 3);
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

	/**
	 * "Position outside, momentum between" variant
	 */
	public void solve () {
		integrator.solve(this);
	}
	
	/**
	 * Writes out current system state as JSON-formatted text
	 * @return the JSON string
	 */
	public String particlesJson () {
		StringBuilder json = new StringBuilder("[");
		Iterator<Particle> p = particles.iterator();
		json.append(p.next().toString());
		while (p.hasNext()) {
			json.append(',' + p.next().toString());
		}
		return (json + "]");
	}
	
	/**
	 * Read initial conditions from a JSON-formatted file
	 * @param fileName the path to the file
	 * @return a Symplectic instance
	 */
	@SuppressWarnings("unchecked")
	public static Symplectic icJson (String fileName) throws IOException {
		BufferedReader bufferedReader = new BufferedReader(new FileReader(new File(fileName)));
		String data = "";
		String line = bufferedReader.readLine();
		while (line != null) {
			data += line;
			line = bufferedReader.readLine();
		}
		bufferedReader.close();
		JSONObject ic = (JSONObject)JSONValue.parse(data);
		List<Particle> bodies = new ArrayList<Particle>();
		for (JSONObject p : (List<JSONObject>)ic.get("bodies")) {
			bodies.add(new Particle((Double)p.get("qX"), (Double)p.get("qY"), (Double)p.get("qZ"), (Double)p.get("pX"), (Double)p.get("pY"), (Double)p.get("pZ"), (Double)p.get("mass")));
		}
		return new Symplectic((Double)ic.get("g"), (Double)ic.get("simulationTime"), (Double)ic.get("timeStep"), (Double)ic.get("errorLimit"), bodies, ((Long)ic.get("integratorOrder")).intValue());
	}
	
	/**
	 * Test method for symplectic integrators
	 * 
	 * @param args None defined
	 * @throws IOException 
	 * @throws Exception 
	 */
	public static void main (String[] args) throws IOException {
		Symplectic scenario;
		if (args.length == 1) {
			scenario = Symplectic.icJson(args[0]);
		} else {
			System.err.println("Missing file name, giving up!");
			return;
		}
		double h0 = scenario.hamiltonian();
		double hMin = h0;
		double hMax = h0;
		System.out.println(scenario.particlesJson());
		System.err.printf("{\"t\":%.2f, \"H\":%.9e, \"H0\":%.9e, \"H-\":%.9e, \"H+\":%.9e, \"ER\":%.1f}%n", 0.0, h0, h0, h0, h0, -999.9);
		long n = 1;
		while (n <= scenario.iterations) {
			scenario.solve();
			double hNow = scenario.hamiltonian();
			double tmp = Math.abs(hNow - h0);
			double dH = tmp > 0.0 ? tmp : 1.0e-18;
			if (hNow < hMin) {
				hMin = hNow;
			} else if (hNow > hMax) {
				hMax = hNow;
			}
			double dbValue = 10.0 * Math.log10(Math.abs(dH / h0));
			System.out.println(scenario.particlesJson());
			System.err.printf("{\"t\":%.2f, \"H\":%.9e, \"H0\":%.9e, \"H-\":%.9e, \"H+\":%.9e, \"ER\":%.1f}%n", n * scenario.timeStep, hNow, h0, hMin, hMax, dbValue);
			if (dbValue > scenario.errorLimit) {
				return;
			}
			n += 1;
		}
	}
}
