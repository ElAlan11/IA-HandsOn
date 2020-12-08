package assignments;

import jade.core.AID;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

public class RegressionGUI extends JFrame {
	//private Agents ag = new Agents();
	private JTextField xField;
	
	
	//Agents a
	RegressionGUI(){
		//super(a.getLocalName());
		//ag = a;
		
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
					String price = xField.getText().trim();
					//ag.updateCatalogue(xValue, Integer.parseInt(price));
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
		addWindowListener(new	WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				//ag.doDelete();
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
