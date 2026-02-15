from koral_input import KORALInput
from koral_settings import KORALSettings

def valid(input_params: KORALInput, settings: KORALSettings) -> bool:
    return (validateInput(input_params) and validateSettings(settings))

def validateInput(input_params: KORALInput) -> bool:
    return True

def validateSettings(settings: KORALSettings) -> bool:
    return True