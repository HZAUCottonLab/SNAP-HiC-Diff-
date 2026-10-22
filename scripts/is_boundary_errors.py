class ISBoundaryError(RuntimeError):
    def __init__(self, code, **context):
        self.code = code
        self.context = context
        details = " ".join(
            f"{key}={value}" for key, value in context.items() if value is not None
        )
        super().__init__(f"ERROR[{code}] {details}".rstrip())
