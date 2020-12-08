import java.util.Random;

public class GeneticAlgorithm {

	public static void main(String[] args) {
		MaximizeFunction mf = new MaximizeFunction(8, 0.6, 0.3, 0);
	}
}

class Chromosome{
	private String chromosome;
	private int fitness;
	
	public Chromosome(String ch) {
		chromosome = ch;
		fitness = calculateFitness();
	}
	
	public void setChromosome(String ch) {
		chromosome = ch;
	}
	
	public String getChromosome(){
		return chromosome;
	}
	
	public int getFitness() {
		return fitness;
	}
	
	private int calculateFitness(){
		int num = Integer.parseUnsignedInt(chromosome, 2);
		int fx = (int) Math.pow(num, 2);
		return fx;
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
			String chromosome = "";
			for(int j=0; j<5; j++) {
				chromosome += Integer.toString(r.nextInt(2));
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

class MaximizeFunction{
	private int populationSize;
	private double mutationRate;
	private double crossoverRate;
	private int elitism;
	private double probabilities[][];
	private Chromosome globalMax;
	private int genOfMax;
	
	private Population population;
	
	public MaximizeFunction(int populSize, double mutRate, double xoverRate, int elitism) {
		this.populationSize = populSize;
		this.mutationRate = mutRate;
		this.crossoverRate = xoverRate;
		this.globalMax = new Chromosome("0");
		
		
		this.elitism = elitism;
		
		population = new Population(populationSize);
		population.initPopulation();
		
		int genCount = 0;

		while(genCount < 10) {
			System.out.println();
			for(Chromosome ch : population.population)
				System.out.println("Chr: " + ch.getChromosome() + " Fitness: " + ch.getFitness());
			
			findMaximum();
			selectParents();
			
			for(double[] p : probabilities)
				System.out.println("Prop: " + p[0] + "\tAcc: " + p[1]);
			
			Population newGeneration = new Population(populationSize);
			crossover(newGeneration);
			mutation(newGeneration);
			
			population = newGeneration;
			genCount++;
		}
		
		System.out.println("\nNumber of generations: " + (genCount) + " generations.");
		System.out.println("Optimal solution = " + globalMax.getChromosome());
	}
	
	public void crossover(Population newGen) {
		Random rGen = new Random();
		
		for(int i=0; i<populationSize; i++) {
			double r = rGen.nextDouble();
			
			if(crossoverRate > r) {
				int crossoverPoint = rGen.nextInt(4) + 1;
				Chromosome secondParent = selectParent();
				String part1 = secondParent.getChromosome().substring(0,crossoverPoint);
				String part2 = population.population[i].getChromosome().substring(crossoverPoint);
				newGen.population[i] = new Chromosome(part1+part2);
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
				int xPoint = rGen.nextInt(5);
				char[] chArray = population.population[i].getChromosome().toCharArray();

				if(chArray[xPoint] == '0') chArray[xPoint] = '1';
				else chArray[xPoint] = '0';
				
				newGen.population[i] = new Chromosome(String.valueOf(chArray));
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
	
	public void findMaximum() {
		for(Chromosome ch : population.population)
			if(ch.getFitness() > globalMax.getFitness())
				globalMax = ch;
	}
}

