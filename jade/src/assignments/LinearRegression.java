/**
 * ***************************************************************
 * JADE - Java Agent DEvelopment Framework is a framework to develop
 * multi-agent systems in compliance with the FIPA specifications.
 * Copyright (C) 2000 CSELT S.p.A.
 * 
 * GNU Lesser General Public License
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA  02111-1307, USA.
 * **************************************************************
 */

package assignments;

import jade.core.Agent;
import jade.core.behaviours.Behaviour;
import jade.core.behaviours.CyclicBehaviour;
import jade.core.behaviours.OneShotBehaviour;
import java.util.Scanner;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

/**
 * This example shows the basic usage of JADE behaviours.<br>
 * More in details this agent executes a <code>CyclicBehaviour</code> that shows
 * a printout at each round and a generic behaviour that performs four successive
 * "dummy" operations. The second operation in particular involves adding a
 * <code>OneShotBehaviour</code>. When the generic behaviour completes the
 * agent terminates.
 * @author Giovanni Caire - TILAB
 */
public class LinearRegression extends Agent {

  private int x; 

  protected void setup() {
    System.out.println("Agent "+getLocalName()+" started.");

    RegressionGUI gui = new RegressionGUI(this);
    gui.showGui();

  } 

  public void updateX(final int xV) {
    addBehaviour(new OneShotBehaviour() {
      public void action() {
        x = xV;
        //addBehaviour(new benettonPred(x));
        addBehaviour(new gradientDescent(x));
      }
    } );
  }

  class benettonPred extends OneShotBehaviour {
    private int n=0;

    benettonPred(int value){
      n = value;
    }
   
    public void action() {
      System.out.println("In Benetton");


      int[][] dataset = {{1,2,3,4,5,6,7,8,9,10},{1,2,3,4,5,6,7,8,9,10},{1,2,3,4,5,6,7,8,9,10}};
      double b0=0, b1=0, yHat=0;
      //int n=0;
      
      //Scanner sc = new Scanner(System.in);
      System.out.println("X = " + n);
      //n = sc.nextInt();
      
      double sumY=0, sumX2=0, sumX=0, sumXY=0;
      
      for(int y : dataset[1])
        sumY += y;
      
      for(int x : dataset[2])
        sumX2 += Math.pow(x, 2);
      
      for(int x : dataset[2])
        sumX += x;
      
      for(int i=0; i<dataset[1].length; i++)
        sumXY += dataset[2][i] * dataset[1][i];
      
      b0 = (sumY*sumX2-sumX*sumXY)/(n*sumX2-Math.pow(sumX, 2));
      b1 = (n*sumXY-sumX*sumY)/(n*sumX2-Math.pow(sumX, 2));
      yHat = b0 + b1*n;
      
      System.out.println("y = b0+b1(x)");
      System.out.println("y = " + b0 + "+" + b1 + "(" + n + ")");
      System.out.println("y = " + yHat);

    }
  }

  class gradientDescent extends OneShotBehaviour {

    private int x;

    gradientDescent(int value){
      x = value;
    } 
   
    public void action() {
      System.out.println("\nIn gradientDescent");

      int[][] dataset = {{36,28,35,39,30,30,31,38,36,38,29,26},{31,29,34,35,29,30,30,38,34,33,29,26}};
      
      int nEpochs = 20;
      double alpha = 0.001;
      double yHat, p, error, b0=0.0, b1=0.0;
      
      //Scanner sc = new Scanner(System.in);
      System.out.println("X = " + x);
      //x = sc.nextInt();
      
      for(int j=0; j<nEpochs; j++) {
        for(int i=0; i<dataset[0].length; i++) {
          p = b0 + b1*dataset[0][i];
          error = p-dataset[1][i];
          
          b0 = b0 - alpha*error;
          b1 = b1 - alpha*error*dataset[0][i];
        }
      }
    
      yHat = b0 + b1*x;
    
      System.out.println("y = b0+b1(x)");
      System.out.println("y = " + b0 + "+" + b1 + "(" + x + ")");
      System.out.println("y = " + yHat);
      

      doDelete();
    }
  }

  class RegressionGUI extends JFrame {
    private LinearRegression ag = new LinearRegression();
    private JTextField xField;
    

    RegressionGUI(LinearRegression a){
      super(a.getLocalName());
      ag = a;
      
      JPanel p = new JPanel();
      p.setLayout(new GridLayout(2, 2));
      p.add(new JLabel("X Value:"));
      xField = new JTextField(15);
      p.add(xField);
      getContentPane().add(p, BorderLayout.CENTER);
      
      JButton addButton = new JButton("Submit");
      addButton.addActionListener( new ActionListener() {
        public void actionPerformed(ActionEvent ev) {
          try {
            String x = xField.getText().trim();
            ag.updateX(Integer.parseInt(x));
            xField.setText("");
          }
          catch (Exception e) {
            JOptionPane.showMessageDialog(RegressionGUI.this, "Invalid values. "+e.getMessage(), "Error", JOptionPane.ERROR_MESSAGE); 
          }
        }
      } );
      p = new JPanel();
      p.add(addButton);
      getContentPane().add(p, BorderLayout.SOUTH);
      
      // Make the agent terminate when the user closes 
      // the GUI using the button on the upper right corner 
      addWindowListener(new WindowAdapter() {
        public void windowClosing(WindowEvent e) {
          ag.doDelete();
          //ag.super.onEnd();
        }
      } );
      
      setResizable(false);
    }
    
    public void showGui() {
      pack();
      Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
      int centerX = (int)screenSize.getWidth() / 2;
      int centerY = (int)screenSize.getHeight() / 2;
      setLocation(centerX - getWidth() / 2, centerY - getHeight() / 2);
      super.setVisible(true);
    }
  } 

}
