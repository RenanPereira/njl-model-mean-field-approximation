# Docker Usage

## Build and Run Locally

Build the Docker image:

Commands to build and locally run image:
```bash
docker build -t njl-model-mean-field-approximation:1.0 .
```

Run the container interactively with the current project directory mounted into the container:
```bash
docker run --rm -it -v $(pwd):/workdir njl-model-mean-field-approximation:1.0 /bin/bash
```

## GitLab Container Registry

Authenticate with GitLab Registry
```bash
docker login registry.gitlab.com
```

Build and Tag the Image
```bash
docker build -t registry.gitlab.com/nambu-jona-lasinio-model/njl-model-mean-field-approximation:1.0 .
```

Push the Image to GitLab Registry
```bash
docker push registry.gitlab.com/nambu-jona-lasinio-model/njl-model-mean-field-approximation:1.0
```

Run the Image from GitLab Registry
```bash
docker run --rm -it -v $(pwd):/workdir registry.gitlab.com/nambu-jona-lasinio-model/njl-model-mean-field-approximation:1.0 /bin/
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