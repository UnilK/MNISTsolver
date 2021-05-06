
This is a personal "learning machine learning the hard way" -project.
Things are designed for maximum experimentation potential.



Here are my initial thoughts on the project:

Words like artif##### ########ence and mach### #####ing
make me sick. My data churning contraptions are called
cakes.

Generally, what a cake does:
It takes in some data, churns it in multiple layers and
then checks the results. Based on how wrong the results
are, the cake is adjusted to do better next time.

Example: feed the cake an image of a 6 and expect it
to classify it as a 6. If it doesn't, hit it with a hammer
and see if it gets the next number right.

Needless to say that this is quite a monkey-banging-rocks-together-ish
approach to solve any problem - but hey, smashing stones
is a fine and fulfilling hobby I would recommend to any primate.

The "how" of cakes is described in the layers that
churn (mush, mix, whatever) the data. A layer is something
that takes in a vector of data, mixes it somehow, and
then passes it on to the next layer. The methods
the layers use are, again, mostly arbitary.

Example: a layer takes in a vector of 256 floating point numbers,
connects to a another layer of size 128 and uses the following
function to do so: layer128[i/2] = layer256[i]-layer256[i+1]

reverse engineerable layers are a good thing to have
for "check results -> what went wrong? -> adjust",
though "check resluts -> favor things that worked
& do some quessing -> adjust" works aswell.

There is no template class for a cake.

I hope I have embraced the vagueness of this all correctly.
