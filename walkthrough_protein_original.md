---
layout: page
title: Protein Walkthrough (original version)
hide_title: true
hide_main_links: true
permalink: /walkthrough_protein_original/
order: 4
---

This walkthrough provides screenshots for procedures described on the main page: [creating simple simulations](/index#run-a-small-simulation) and [running a simple analysis](/index#run-a-simple-analysis).

Welcome to the simulator.

![alt](../images/walkthrough/w1-simulator.png)

Make a new simulation. Choose a relevant name.

![alt](../images/walkthrough/w2-simulation-new.png)

Once you make a new simulation you choose a kickstarter. For this example select the "proteins".

![alt](../images/walkthrough/w3-kickstart.png)

Choose an experiment. Select "protein".

![alt](../images/walkthrough/w4-experiment.png)

Use the default settings. And start the simulation.

![alt](../images/walkthrough/w5-settings.png)

Back on the simulator page, we can see which simulations are running.

![alt](../images/walkthrough/w6-running.png)

On the simulation page, you can watch the AUTOMACS log.

![alt](../images/walkthrough/w7-logging.png)

Back on the calculator page. Make the default meta file.

![alt](../images/walkthrough/w8-calculator.png)

You can also inspect the simulation times.

![alt](../images/walkthrough/w9-generate-meta.png)

The "inspect times" notebook prints out your simulations and the start and end times for each part.

![alt](../images/walkthrough/w9b-inspect-times.png)

You can also view the automatically generated meta file, which is always called "meta.current.yaml" by clicking on the link on the main page.

![alt](../images/walkthrough/w9c-inspect-times-notebook.png)

To run a computation, select the desired meta file and then hit "compute".

![alt](../images/walkthrough/w10-view-meta.png)

The calculation is logged to a console. First we make simulation slices and then we run the "protein_rmsd" calculation.

![alt](../images/walkthrough/w11-compute.png)

Once the calculation is complete, click on "protein_rmsd" in the plots tile.

![alt](../images/walkthrough/w12-computing.png)

You can generate an interactive notebook from the plot script. 

![alt](../images/walkthrough/w13-make-notebook.png)

If you run the notebook it creates a plot which is also saved to the omnicalc plot location.

![alt](../images/walkthrough/w14-notebook-plot.png)

Back on the calculator page, you should "regenerate the thumbnails" and refresh the page.

![alt](../images/walkthrough/w15-regenerate-thumbnails.png)

Then you can show the pictures using the toggle buttons.

![alt](../images/walkthrough/w16-view-pictures.png)

After you run longer trajectories from the terminal (or on a cluster), you can regenerate the meta file to extend the trajectory and analyze a longer simulation.

![alt](../images/walkthrough/w17-remake-meta.png)

More data means more plots.

![alt](../images/walkthrough/w18-wow.png)

