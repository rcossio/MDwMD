const express = require('express');
const bodyParser = require('body-parser');
const mongoose = require('mongoose');
const jwt = require('jsonwebtoken');
const cookieParser = require('cookie-parser');
const config = require('./config/environment');

const mainRouter = require('./routes/main');
const authRouter = require('./routes/auth');

const app = express();
const PORT = config.port;

// Connect to MongoDB
mongoose.connect(config.mongo.uri)
.then(() => console.log('Connected to MongoDB'))
.catch(err => console.error('Error connecting to MongoDB:', err));

// Middleware
app.use(bodyParser.json());
app.use(cookieParser());
app.use(express.static('public'));

// Middleware to verify JWT token
const verifyToken = (req, res, next) => {
  const token = req.cookies.jwt || ''; // Assuming the token is sent in a cookie
  jwt.verify(token, config.jwtSecret, (err, decoded) => {
    if (err) {
      return res.status(401).json({ message: 'Unauthorized access' });
    }
    req.user = decoded; // Add decoded token to request so it can be used in the route
    next();
  });
};

// A protected route
app.get('/protected', verifyToken, (req, res) => {
  res.json({ message: 'This is a protected route', user: req.user });
});

// Use router
app.use('/', mainRouter);
app.use('/auth', authRouter);


app.listen(PORT, () => {
  console.log(`Server running on http://localhost:${PORT}`);
});
