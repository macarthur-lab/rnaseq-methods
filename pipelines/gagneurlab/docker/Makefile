TAG1 = weisburd/gagneurlab:v1.1.7
TAG2 = weisburd/gagneurlab:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
