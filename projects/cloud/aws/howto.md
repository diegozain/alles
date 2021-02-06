# How to AWS

Taken from [here](https://machinelearningmastery.com/develop-evaluate-large-deep-learning-models-keras-amazon-web-services/) and [here](https://machinelearningmastery.com/command-line-recipes-deep-learning-amazon-web-services/).

Also look at the official manpage on [EC2](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/concepts.html).

1. go to https://console.aws.amazon.com/console/home

1. choose ```EC2 Virtual Servers in the Cloud ```.

1. Top right choose ```US West Oregon```.

1. Click ```Launch instance```.

1. Go to ```Community AMIs``` and in the search bar look for ```deep learning ami```.

1. Choose ```Deep Learning AMI (Amazon Linux)```.

1. Choose ```p3.2xlarge```.

1. Click ```Review and Launch``` and then ```Launch```. Make sure you have the key-pair.

1. Click ```View Instances```.

1. __You're set now__. Look for the public IP in the ```View Instances``` dashboard.

1. Open a terminal where the key-pair is and type:

	```ssh -i keras-keypair.pem ec2-user@54.186.97.77```

1. You now have to choose your environment bs:

	```source activate tensorflow_p36```

1. Do whatever you want, and then type 

	```exit```

1. Go back to the AWS console webpage, click on ```EC2``` and choose ```Instances```.

1. Click the ```Actions``` button, then ```Instance State```, then ```Terminate```.

### SCP to AWS

```scp -i keras-keypair.pem file.py ec2-user@54.186.97.77:~/```

### SCP from AWS

```scp -i keras-keypair.pem ec2-user@54.186.97.77:~/file.png```

### Run job in AWS and close local terminal

```nohup python script.py >output.log </dev/null 2>&1 &```

### Stop AWS from closing inactive local terminal

```watch "tail output.log"```

### Use discarded space - Spot Instances

1. These instances are subject to unannounced termination at any moment!  

1. Go to ```Spot Requests``` under ```Instances```.

1. Jobs should be:

	* Instance type flexible
	* Interruption tolerant
	* Use dynamic scaling prices

1. s
