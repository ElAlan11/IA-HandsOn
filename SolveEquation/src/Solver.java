import java.util.Random;
import java.util.Arrays;

public class Solver {
	
	public static void main(String[] args) {
		GeneticAlgorithm ga = new GeneticAlgorithm(8, 0.6, 0.3, 0);
	}
}

// a + 2b - 3c + d + 4e + f = 30

class Chromosome{
	private int[] chromosome;
	private int fitness;
	
	public Chromosome(int[] ch) {
		chromosome = ch;
		fitness = calculateFitness();
	}
	
	public void setChromosome(int[] ch) {
		chromosome = ch;
		calculateFitness();
	}
	
	public int[] getChromosome(){
		return chromosome;
	}
	
	public int getFitness() {
		return fitness;
	}
	
	private int calculateFitness(){
		int fitness = 100 - Math.abs((chromosome[0] + 2*chromosome[1] - 3*chromosome[2] + chromosome[3] + 4*chromosome[4] + chromosome[5]) - 31);
		return fitness;
	}
}

class Population{
	public Chromosome[] population;
	private int populationFitness;
	private int size;
	
	public Population(int size) {
		population = new Chromosome[size];
		this.size = size;
	}
	
	public int getPopulationFitness() {
		return populationFitness;
	}
	
	public void initPopulation() {
		Random r = new Random();
		
		for(int i=0; i<size; i++) {
			int[] chromosome = new int[6];
			for(int j=0; j<6; j++) {
				chromosome[j] = r.nextInt(10);
			}
			population[i] = new Chromosome(chromosome);
		}
	}
	
	public void calculatePopulationFitness() {
		int sum = 0;

		for(Chromosome ch : population)
			sum += ch.getFitness();
		populationFitness = sum;
	}
}

class GeneticAlgorithm{
	private int populationSize;
	private double mutationRate;
	private double crossoverRate;
	private int elitism;
	private double probabilities[][];
	
	private Population population;
	
	public GeneticAlgorithm(int populSize, double mutRate, double xoverRate, int elitism) {
		this.populationSize = populSize;
		this.mutationRate = mutRate;
		this.crossoverRate = xoverRate;
		
		this.elitism = elitism;
		
		population = new Population(populationSize);
		population.initPopulation();
		
		int genCount = 0;

		while(genCount < 30) {
			System.out.println();
			for(Chromosome ch : population.population)
				System.out.println("Chr: " + Arrays.toString(ch.getChromosome()) + " Fitness: " + ch.getFitness());
			
			Chromosome res = evaluatePopulation();
			if(res != null) {
				System.out.println("\nSolution found after " + (genCount+1) + " generations.");
				System.out.println("Optimal solution = " + Arrays.toString(res.getChromosome()));
				break;
			}

			selectParents();
			
			for(double[] p : probabilities)
				System.out.println("Prop: " + p[0] + "\tAcc: " + p[1]);
			
			Population newGeneration = new Population(populationSize);
			crossover(newGeneration);
			mutation(newGeneration);
			
			population = newGeneration;
			genCount++;
		}
	}
	
	public void crossover(Population newGen) {
		Random rGen = new Random();
		
		for(int i=0; i<populationSize; i++) {
			double r = rGen.nextDouble();
			
			if(crossoverRate > r) {
				int crossoverPoint = rGen.nextInt(5) + 1;
				Chromosome secondParent = selectParent();
				
				int[] res = new int[6];
				System.arraycopy( secondParent.getChromosome(), 0, res, 0, crossoverPoint );
				System.arraycopy( population.population[i].getChromosome(), crossoverPoint, res, crossoverPoint, 6-crossoverPoint );
				newGen.population[i] = new Chromosome(res);
			}
			else
				newGen.population[i] = new Chromosome(population.population[i].getChromosome());
		}
	}
	
	public void mutation(Population newGen) {
		Random rGen = new Random();
		
		for(int i=0; i<populationSize; i++) {
			double r = rGen.nextDouble();
			
			if(mutationRate > r) {
				int xPoint = rGen.nextInt(6);
				int rNum = rGen.nextInt(10);
				int[] newChromosome = new int[6];
				
				System.arraycopy(population.population[i].getChromosome(), 0, newChromosome, 0, 6);
				newChromosome[xPoint] = rNum;
				
				newGen.population[i] = new Chromosome(newChromosome);
			}
			else
				newGen.population[i] = new Chromosome(population.population[i].getChromosome());
		}
	}
	
	public Chromosome selectParent() {
		Random rGen = new Random();
		Chromosome selected;
		double r = rGen.nextDouble();
		int i=0;

		for(; i<populationSize; i++) {
			if(r < probabilities[i][1])
				break;
		}
		
		selected = population.population[i];
		return selected;
	}
	
	public void selectParents() {
		Random rGen = new Random();
		probabilities = new double[populationSize][2];
		double acc = 0;
		Population selectedParents = new Population(populationSize);
		int chromCount = 0;
		
		population.calculatePopulationFitness();
		System.out.println("Population Fitness: " + population.getPopulationFitness());
		
		for(int i=0; i<populationSize; i++) {
			probabilities[i][0] = (double)population.population[i].getFitness() / population.getPopulationFitness();
			acc += probabilities[i][0];
			probabilities[i][1] = acc;
		}
		
		while(chromCount < populationSize) {
			double r = rGen.nextDouble();
			int i=0;
			
			for(; i<populationSize; i++) {
				if(r < probabilities[i][1])
					break;
			}
			
			selectedParents.population[chromCount] = population.population[i];
			chromCount++;
		}
		
		population = selectedParents;
	}
	
	public Chromosome evaluatePopulation() {
		for(Chromosome ch : population.population)
			if(ch.getFitness() == 100)
				return ch;
		return null;
	}
}


