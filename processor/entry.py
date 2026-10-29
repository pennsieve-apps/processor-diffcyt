"""
One image, two runtimes: as a Lambda function (AWS_LAMBDA_RUNTIME_API is set)
it hands over to the Lambda runtime interface client with processor.handler;
as a container it replaces itself with `Rscript R/main.R`, which reads
INPUT_DIR, OUTPUT_DIR and the parameters from the environment.
"""
import os
import sys

from processor.handler import main_script

if os.environ.get("AWS_LAMBDA_RUNTIME_API"):
    os.execvp(sys.executable, [sys.executable, "-m", "awslambdaric", "processor.handler.handler"])

os.execvp("Rscript", ["Rscript", main_script()])
