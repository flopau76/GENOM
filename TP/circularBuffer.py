class CircularBuffer:
    def __init__(self, elements):
        """
        Initialize a circular buffer with a fixed capacity.
        
        :param capacity: Maximum number of elements the buffer can hold
        :param elements: Optional initial list of elements to populate the buffer
        """
        elements = list(elements)
        if len(elements)<= 0:
            raise ValueError("Capacity must be a positive integer")
        
        self.buffer = elements[:]
        self.head = 0

    def enqueue(self, item):
        """
        Add an item to the buffer.      
        :param item: Element to be added
        """
        self.buffer[self.head] = item
        self.head = (self.head +1)%len(self)

    def peek(self):
        """
        View the oldest item in the buffer without removing it.
        
        :return: The oldest item in the buffer
        """
        return self.buffer[self.head]

    def __len__(self):
        """
        Return the current number of elements in the buffer.
        
        :return: Current size of the buffer
        """
        return len(self.buffer)

    def __str__(self):
        """
        Return a string representation of the buffer's contents.
        
        :return: String representation of buffer contents
        """
        
        contents = []
        current = self.head
        for _ in range(self.size):
            contents.append(str(self.buffer[current]))
            current = (current + 1) % len(self)
        
        return f"[{', '.join(contents)}]"
