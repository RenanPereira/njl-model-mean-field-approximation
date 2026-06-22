# Docker Usage

## Build and Run Locally

Build the Docker image:

Commands to build the image:
```bash
docker build -f Dockerfile -t njl-model-mean-field-approximation:1.0 .

docker build -f Dockerfile.gsl-source -t njl-model-mean-field-approximation:1.0-gsl-2.8 .
```
On the first Dockerfile gsl is installed via the system package (`libgsl-dev`), while the second is built from the source code itself stored inside the file `third_party/gsl-2.8.tar.gz`. 


Run the container interactively with the current project directory mounted into the container:
```bash
docker run --rm -it -v $(pwd):/workdir njl-model-mean-field-approximation:1.0 /bin/bash

docker run --rm -it -v $(pwd):/workdir njl-model-mean-field-approximation:1.0-gsl-2.8 /bin/bash
```

## GitLab Container Registry

Authenticate with GitLab Registry
```bash
docker login registry.gitlab.com
```

Build and Tag the Image
```bash
docker build -f Dockerfile -t registry.gitlab.com/nambu-jona-lasinio-model/njl-model-mean-field-approximation:1.0 .

docker build -f Dockerfile.gsl-source -t registry.gitlab.com/nambu-jona-lasinio-model/njl-model-mean-field-approximation:1.0-gsl-2.8 .
```

Push the Image to GitLab Registry
```bash
docker push registry.gitlab.com/nambu-jona-lasinio-model/njl-model-mean-field-approximation:1.0

docker push registry.gitlab.com/nambu-jona-lasinio-model/njl-model-mean-field-approximation:1.0-gsl-2.8
```

Run the Image from GitLab Registry
```bash
docker run --rm -it -v $(pwd):/workdir registry.gitlab.com/nambu-jona-lasinio-model/njl-model-mean-field-approximation:1.0 /bin/bash

docker run --rm -it -v $(pwd):/workdir registry.gitlab.com/nambu-jona-lasinio-model/njl-model-mean-field-approximation:1.0-gsl-2.8 /bin/bash
```

Notes
The `--rm` flag automatically removes the container after exit.
The `-it` flags enable interactive terminal usage.
The volume mount:
```bash
-v $(pwd):/workdir
```

maps the current local repository into the container at `/workdir`.

Any files created inside `/workdir` in the container will persist on the host filesystem.
